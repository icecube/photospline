#!/usr/bin/env python

import pickle
import unittest
import photospline

import numpy as np

from pathlib import Path


class TestEvaluation(unittest.TestCase):
    def setUp(self):
        self.testdata = Path(__file__).parent / "test_data"
        self.spline = photospline.SplineTable(self.testdata / "test_spline_4d.fits")
        extents = np.array(self.spline.extents)
        loc = extents[:, :1]
        scale = np.diff(extents, axis=1)
        self.x = (np.random.uniform(0, 1, size=(self.spline.ndim, 10)) + loc) * scale

    def test_pickle(self):
        restored = pickle.loads(pickle.dumps(self.spline))
        np.testing.assert_allclose(self.spline.evaluate_simple(self.x), restored.evaluate_simple(self.x), rtol=0, atol=0)

if __name__ == "__main__":
    unittest.main()
