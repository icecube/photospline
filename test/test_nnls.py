import unittest
import numpy as np
import photospline


class TestNNLS(unittest.TestCase):
    def setUp(self):
        try:
            from scipy import sparse
        except ImportError:
            raise unittest.SkipTest("test requires scipy")

        self.A = np.array([[1, 0], [1, 0], [0, 1]])
        self.Asp = sparse.csc_matrix(self.A)
        self.b = np.array([2, 1, 1])

    def testSparseIsSparse(self):
        self.assertEqual(len(self.Asp.data), 3)

    def testImplementationMatchesScipy(self):
        from scipy import optimize

        np.testing.assert_allclose(
            photospline.nnls(self.Asp, self.b),
            optimize.nnls(self.A, self.b)[0],
            err_msg="sparse result does not aggree with dense",
        )


if __name__ == "__main__":
    unittest.main()
