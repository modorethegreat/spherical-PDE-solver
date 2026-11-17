"""An exact reference and timestep refinement study for one real harmonic."""

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class ModeCase:
    degree: int = 4
    order: int = 2
    steps: int = 100
    drag: float = 0.5
    duration: float = 1.0
    mean_depth: float = 2.0
    gravity: float = 1.0
    radius: float = 1.0

    def __post_init__(self):
        if not (0 <= self.order <= self.degree <= 100):
            raise ValueError("Require 0 ≤ order ≤ degree ≤ 100.")
        if not 1 <= self.steps <= 20_000:
            raise ValueError("Steps must be between 1 and 20,000.")
        if not (math.isfinite(self.drag) and 0 <= self.drag <= 1e6):
            raise ValueError("Drag must be finite and between 0 and 1,000,000.")
        for name in ("duration", "mean_depth", "gravity", "radius"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive.")
        a, c = self.coupling
        if not all(map(math.isfinite, (a, c, a * c))) or c == 0:
            raise ValueError("The physical parameters exceed the supported floating-point range.")
        if self.duration / self.steps * math.sqrt(a * c) >= 2:
            raise ValueError("Choose more steps so that Δt × wave frequency is below 2.")

    @property
    def coupling(self):
        a = (self.mean_depth / self.radius) * self.degree * (self.degree + 1)
        return a, self.gravity / self.radius


def numerical_step(h, v, a, c, drag, dt):
    """Advance the velocity first with implicit drag, then update height."""
    v_next = (v - dt * c * h) / (1 + dt * drag)
    return h + dt * a * v_next, v_next


def exact_state(t, a, c, drag, h0=1.0, v0=0.0):
    """Apply the exact 2×2 matrix exponential in all damping regimes."""
    half_drag = drag / 2
    discriminant = half_drag**2 - a * c
    if discriminant < 0:
        frequency = math.sqrt(-discriminant)
        decay = math.exp(-half_drag * t)
        diagonal = decay * math.cos(frequency * t)
        off_diagonal = decay * math.sin(frequency * t) / frequency
    elif discriminant == 0:
        diagonal = math.exp(-half_drag * t)
        off_diagonal = t * diagonal
    else:
        root = math.sqrt(discriminant)
        # Re-express the slow eigenvalue to avoid subtracting similar numbers.
        slow = math.exp(-(a * c / (half_drag + root)) * t)
        gap = -math.expm1(-2 * root * t)
        diagonal = slow * (1 - gap / 2)
        off_diagonal = slow * gap / (2 * root)
    return (
        diagonal * h0 + off_diagonal * (half_drag * h0 + a * v0),
        diagonal * v0 - off_diagonal * (c * h0 + half_drag * v0),
    )


def benchmark(case):
    """Return synchronized coefficients and a four-level final-time study."""
    a, c = case.coupling
    weight = math.sqrt(a / c)
    dt = case.duration / case.steps

    def state_error(h, v, h_exact, v_exact):
        return math.hypot(h - h_exact, weight * (v - v_exact))

    h, v = 1.0, 0.0
    samples = []
    for index in range(case.steps + 1):
        t = case.duration if index == case.steps else index * dt
        he, ve = exact_state(t, a, c, case.drag)
        samples.append([t, h, he, v, ve, state_error(h, v, he, ve)])
        if index < case.steps:
            h, v = numerical_step(h, v, a, c, case.drag, dt)

    rows = []
    previous = None
    for level in range(4):
        steps = case.steps * 2**level
        step_size = case.duration / steps
        if level == 0:
            error = samples[-1][-1]
        else:
            h, v = 1.0, 0.0
            for _ in range(steps):
                h, v = numerical_step(h, v, a, c, case.drag, step_size)
            error = state_error(h, v, he, ve)
        observed = math.log2(previous / error) if previous and error else None
        rows.append({"steps": steps, "dt": step_size, "error": error, "order": observed})
        previous = error

    return {
        "schema_version": 1,
        "case": asdict(case),
        "initial_state": {"height": 1.0, "velocity": 0.0},
        "wave_frequency": math.sqrt(a * c),
        "sample_columns": ["time", "height_numerical", "height_exact", "velocity_numerical", "velocity_exact", "state_error"],
        "samples": samples,
        "convergence": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=4)
    parser.add_argument("--order", type=int, help="Default: the lesser of 2 and the degree.")
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--drag", type=float, default=0.5)
    parser.add_argument("--duration", type=float, default=1.0)
    parser.add_argument("--mean-depth", type=float, default=2.0)
    parser.add_argument("--gravity", type=float, default=1.0)
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--output", type=Path, help="Write JSON to a file instead of stdout.")
    args = parser.parse_args()
    try:
        result = benchmark(ModeCase(
            args.degree, min(2, args.degree) if args.order is None else args.order,
            args.steps, args.drag, args.duration, args.mean_depth, args.gravity, args.radius,
        ))
        payload = json.dumps(result, indent=2, allow_nan=False)
    except ValueError as error:
        parser.error(str(error))
    if args.output:
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)


if __name__ == "__main__":
    main()
