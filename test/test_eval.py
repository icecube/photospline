#!/usr/bin/env python

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

    def test_vector(self):
        centers = self.spline.search_centers(self.x)
        self.assertIsInstance(centers, np.ndarray)
        self.assertEqual(centers.shape, (self.spline.ndim, self.x.shape[1]))
        v = self.spline.evaluate(self.x, centers=centers)
        self.assertIsInstance(v, np.ndarray)
        self.assertEqual(v.shape, (self.x.shape[1],))

    def test_scalar(self):
        centers = self.spline.search_centers([x[0] for x in self.x])
        self.assertIsInstance(centers, tuple)
        self.assertEqual(len(centers), self.spline.ndim)
        v = self.spline.evaluate([x[0] for x in self.x], centers=centers)
        self.assertIsInstance(v, float)


if __name__ == "__main__":
    unittest.main()
