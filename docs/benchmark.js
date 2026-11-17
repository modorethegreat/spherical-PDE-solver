// This browser port is checked against mode_benchmark.py in tests/test_browser_parity.py.
export function modeBenchmark(degree, order, steps, drag) {
  if (!Number.isInteger(degree) || degree < 1 || degree > 12 ||
      !Number.isInteger(order) || order < 0 || order > degree ||
      !Number.isInteger(steps) || steps < 20 || steps > 640 ||
      !Number.isFinite(drag) || drag < 0 || drag > 4) {
    throw new RangeError('Choose a valid degree, order, timestep count and drag coefficient.');
  }

  const duration = 1, meanDepth = 2, gravity = 1, radius = 1;
  const coupling = meanDepth / radius * degree * (degree + 1);
  const pressure = gravity / radius;
  const frequency = Math.sqrt(coupling * pressure);
  const velocityWeight = Math.sqrt(coupling / pressure);
  const halfDrag = drag / 2;
  const discriminant = halfDrag ** 2 - frequency ** 2;
  const dt = duration / steps;

  function exact(t) {
    let diagonal, offDiagonal;
    if (discriminant < 0) {
      const omega = Math.sqrt(-discriminant);
      const decay = Math.exp(-halfDrag * t);
      diagonal = decay * Math.cos(omega * t);
      offDiagonal = decay * Math.sin(omega * t) / omega;
    } else if (discriminant === 0) {
      diagonal = Math.exp(-halfDrag * t);
      offDiagonal = t * diagonal;
    } else {
      const root = Math.sqrt(discriminant);
      const slow = Math.exp(-coupling * pressure * t / (halfDrag + root));
      const gap = -Math.expm1(-2 * root * t);
      diagonal = slow * (1 - gap / 2);
      offDiagonal = slow * gap / (2 * root);
    }
    return [diagonal + halfDrag * offDiagonal, -pressure * offDiagonal];
  }

  function trajectory(count, visit) {
    const step = duration / count;
    let h = 1, v = 0;
    visit?.(0, h, v);
    for (let index = 1; index <= count; index++) {
      v = (v - step * pressure * h) / (1 + step * drag);
      h += step * coupling * v;
      visit?.(index === count ? duration : index * step, h, v);
    }
    return [h, v];
  }

  const error = (h, v, he, ve) => Math.hypot(h - he, velocityWeight * (v - ve));
  const samples = [];
  trajectory(steps, (t, h, v) => {
    const [he, ve] = exact(t);
    samples.push([t, h, he, v, ve, error(h, v, he, ve)]);
  });

  const final = exact(duration);
  const convergence = [];
  let previous;
  for (let level = 0; level < 4; level++) {
    const count = steps * 2 ** level;
    const mismatch = level === 0 ? samples.at(-1)[5] : error(...trajectory(count), ...final);
    convergence.push({ steps: count, dt: dt / 2 ** level, error: mismatch,
      order: previous && mismatch ? Math.log2(previous / mismatch) : null });
    previous = mismatch;
  }

  return {
    schema_version: 1,
    case: { degree, order, steps, drag, duration, mean_depth: meanDepth, gravity, radius },
    initial_state: { height: 1, velocity: 0 },
    wave_frequency: frequency,
    sample_columns: ['time', 'height_numerical', 'height_exact', 'velocity_numerical', 'velocity_exact', 'state_error'],
    samples, convergence,
  };
}
