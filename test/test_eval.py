#!/usr/bin/env python

import sys, os

sys.path.append(os.getcwd())
import unittest
import photospline

import numpy as np

from pathlib import Path

class TestEvaluation(unittest.TestCase):

    def setUp(self):
        self.testdata = Path(__file__).parent / "test_data"
        self.spline = photospline.SplineTable(self.testdata / "test_spline_4d.fits")
        extents = np.array(self.spline.extents)
        loc = extents[:,:1]
        scale = np.diff(extents, axis=1)
        self.x = (np.random.uniform(0, 1, size=(self.spline.ndim, 10)) + loc) * scale

    def test_searchcenters_vector(self):
        centers = self.spline.search_centers(self.x)
        self.assertIsInstance(centers, np.ndarray)
        self.assertEqual(centers.shape, (self.spline.ndim, self.x.shape[1]))
    
    def test_searchcenters_scalar(self):
        centers = self.spline.search_centers([x[0] for x in self.x])
        self.assertIsInstance(centers, tuple)
        self.assertEqual(len(centers), self.spline.ndim)


if __name__ == "__main__":
    unittest.main()