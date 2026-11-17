import unittest

import numpy as np

from mode_benchmark import ModeCase, exact_state, numerical_step
from sph_harm_transform import ISHT, SHT, Y, ls, phi, weights


class SpectralTests(unittest.TestCase):
    def test_real_zonal_cosine_and_sine_modes(self):
        for degree, order, quadrature, expected in [
            (3, 0, "real", 1), (4, 2, "real", 0.5), (4, 2, "imag", -0.5j),
        ]:
            with self.subTest(degree=degree, order=order, quadrature=quadrature):
                field = getattr(Y(degree, order), quadrature)
                coefficients = SHT(field)
                expected_coefficients = np.zeros_like(coefficients)
                expected_coefficients[degree, order] = expected
                np.testing.assert_allclose(coefficients, expected_coefficients, atol=8e-13)
                np.testing.assert_allclose(ISHT(coefficients), field, atol=1e-11)

    def test_grid_to_grid_convergence(self):
        case = ModeCase()
        field = Y(case.degree, case.order).real
        h0 = SHT(field)
        coupling = (case.mean_depth / case.radius) * (ls * (ls + 1))[:, None]
        time = 0.4
        expected = exact_state(time, *case.coupling, case.drag)[0] * field
        errors = []
        for steps in (25, 50, 100):
            h, v = h0.copy(), np.zeros_like(h0)
            for _ in range(steps):
                h, v = numerical_step(h, v, coupling, case.gravity / case.radius, case.drag, time / steps)
            error = ISHT(h) - expected
            errors.append(np.sqrt(np.sum(weights[:, None] * error**2) * (2*np.pi / len(phi))))
        for coarse, fine in zip(errors, errors[1:]):
            self.assertGreater(coarse / fine, 1.8)
            self.assertLess(coarse / fine, 2.2)


if __name__ == "__main__":
    unittest.main()
