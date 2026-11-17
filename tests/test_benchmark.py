import math
import unittest

from mode_benchmark import ModeCase, benchmark, exact_state, numerical_step


class ModeBenchmarkTests(unittest.TestCase):
    def test_undamped_reference_from_sine_and_cosine(self):
        case = ModeCase(degree=4, order=2, drag=0)
        a, c = case.coupling
        frequency = math.sqrt(a * c)
        t = 0.27
        initial_height, initial_velocity = 0.7, -0.12
        h, v = exact_state(t, a, c, 0, initial_height, initial_velocity)
        self.assertAlmostEqual(h, initial_height * math.cos(frequency * t)
                               + a / frequency * initial_velocity * math.sin(frequency * t))
        self.assertAlmostEqual(v, initial_velocity * math.cos(frequency * t)
                               - c / frequency * initial_height * math.sin(frequency * t))

    def test_critical_and_overdamped_reference(self):
        a = c = 1.0
        h, v = exact_state(0.5, a, c, 2.0)
        self.assertAlmostEqual(h, 1.5 * math.exp(-0.5))
        self.assertAlmostEqual(v, -0.5 * math.exp(-0.5))
        h, v = exact_state(1.0, a, c, 2.5)
        self.assertAlmostEqual(h, (4 * math.exp(-0.5) - math.exp(-2)) / 3)
        self.assertAlmostEqual(v, (-2 * math.exp(-0.5) + 2 * math.exp(-2)) / 3)

    def test_exact_energy_budget(self):
        a, c = ModeCase().coupling
        for drag in (0, 0.5):
            dt = 0.001
            states = [exact_state(i * dt, a, c, drag) for i in range(1001)]
            initial_energy = 0.5 * c
            h, v = states[-1]
            final_energy = 0.5 * (c * h**2 + a * v**2)
            dissipated = drag * a * dt * (
                sum(v**2 for _, v in states[1:-1]) + 0.5 * states[-1][1]**2
            )
            self.assertAlmostEqual(initial_energy, final_energy + dissipated, places=6)

    def test_the_mean_height_is_unchanged(self):
        for drag in [0, 0.5]:
            self.assertEqual(exact_state(20, 0, 1, drag, 2, 0)[0], 2)
            self.assertEqual(numerical_step(2, -1, 0, 1, drag, 0.1)[0], 2)

    def test_initial_state_and_first_order_refinement(self):
        result = benchmark(ModeCase(degree=4, order=2, steps=100))
        self.assertEqual(result["samples"][0], [0, 1, 1, 0, 0, 0])
        self.assertEqual(result["samples"][-1][0], 1)
        rows = result["convergence"]
        self.assertEqual([row["steps"] for row in rows], [100, 200, 400, 800])
        for row in rows[1:]:
            self.assertGreater(row["order"], 0.85)
            self.assertLess(row["order"], 1.05)

    def test_invalid_parameters_and_unstable_step(self):
        for params in [{"degree": 2, "order": 3}, {"drag": float("nan")},
                       {"steps": 0}, {"degree": 100, "steps": 2}]:
            with self.subTest(params=params), self.assertRaises(ValueError):
                ModeCase(**params)


if __name__ == "__main__":
    unittest.main()
