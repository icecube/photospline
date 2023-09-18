#!/usr/bin/env python

import sys, os

sys.path.append(os.getcwd())
import numpy
import photospline

from pathlib import Path

import unittest


class StackTest(unittest.TestCase):
    def setUp(self):
        self.testdata = Path(__file__).parent / "test_data"
        self.spline = photospline.SplineTable(
            str(self.testdata / "test_spline_1d.fits")
        )

    def testStack(self):
        spline = self.spline
        stack = photospline.SplineTable.stack([spline] * 3, [-1, 0, 1], stackOrder=2)
        self.assertEqual(spline.ndim, 1)
        self.assertEqual(stack.ndim, 2)
        self.assertEqual(stack.order, (2, 2))
        numpy.testing.assert_equal(stack.knots[0], spline.knots[0])
        numpy.testing.assert_equal(
            stack.knots[1], [-3.6, -2.6, -1.6, -0.6, 0.4, 1.4, 2.4, 3.4]
        )

    def testBadArgs(self):
        with self.assertRaises(ValueError):
            photospline.SplineTable.stack([], [])
        with self.assertRaises(TypeError):
            photospline.SplineTable.stack(range(3), [])
        with self.assertRaises(TypeError):
            photospline.SplineTable.stack([self.spline], ["a"])
        with self.assertRaises(TypeError):
            photospline.SplineTable.stack([self.spline, None], [])
        with self.assertRaises(ValueError):
            photospline.SplineTable.stack([self.spline] * 2, [1])

    def testSetup(self):
        assert self.testdata.exists()


if __name__ == "__main__":
    unittest.main()
