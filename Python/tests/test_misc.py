import unittest
import operator
import numpy as np
import numpy.testing
import math

from src.prcv.misc import get_cvncomp, get_cvsettings

class TestCVMethods(unittest.TestCase):

    def test_errors(self):
        """
        Tests for error handlings.
        """

        # wrong type of "cv" argument
        self.assertRaises(ValueError, get_cvsettings, (1), 1, 1)

        # mismatch "cv" between size of "cv" and "nrows"
        self.assertRaises(ValueError, get_cvsettings, np.array([1, 2, 1, 2, 1, 2]), 7, 1)
        self.assertRaises(ValueError, get_cvsettings, np.array([1, 2, 1, 2, 1, 2]), 5, 1)

        # mismatch "cv" between size of "cv" and "nrows"
        self.assertRaises(ValueError, get_cvsettings, {'cv': 'rand', 'nseg': 2}, 5, 1)
        self.assertRaises(ValueError, get_cvsettings, {'type': 'rand', 'nseg': 1}, 5, 1)
        self.assertRaises(ValueError, get_cvsettings, {'type': 'rand', 'nseg': 6}, 5, 1)
        self.assertRaises(ValueError, get_cvsettings, {'type': 'ven', 'nseg': 2}, 6, 1, np.array([1.0, 2.0]))
        self.assertRaises(ValueError, get_cvsettings, {'type': 'ven', 'nseg': 2}, 4, 1, np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
        self.assertRaises(ValueError, get_cvsettings, {'type': 'grand', 'nseg': 2}, 10, 1)


    def test_get_ncomp(self):
        """
        Test that number of components is equal to number of rows in local
        calibration set minus 1.
        """

        # A = 1
        self.assertEqual(get_cvncomp(10, 1, 2), 1)
        self.assertEqual(get_cvncomp(10, 1, 5), 1)
        self.assertEqual(get_cvncomp(10, 1, 6), 1)
        self.assertEqual(get_cvncomp(10, 1, 10), 1)

        # A = 5
        self.assertEqual(get_cvncomp(10, 5, 2), 4)
        self.assertEqual(get_cvncomp(10, 5, 5), 5)
        self.assertEqual(get_cvncomp(10, 5, 6), 5)
        self.assertEqual(get_cvncomp(10, 5, 10), 5)

        # A = 9
        self.assertEqual(get_cvncomp(10, 9, 2), 4)
        self.assertEqual(get_cvncomp(10, 9, 5), 7)
        self.assertEqual(get_cvncomp(10, 9, 6), 7)
        self.assertEqual(get_cvncomp(10, 9, 10), 8)

        # A = 10
        self.assertEqual(get_cvncomp(10, 10, 2), 4)
        self.assertEqual(get_cvncomp(10, 10, 5), 7)
        self.assertEqual(get_cvncomp(10, 10, 6), 7)
        self.assertEqual(get_cvncomp(10, 10, 10), 8)


    def test_get_settings_manual(self):
        """
        Test that "get_settings" return manually defined vector with indices correctly.
        """

        # A = 5, K = 3
        ind = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings(ind, 10, 5)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 7, K = 3
        ind = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings(ind, 10, 7)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 7, K = 4
        ind = np.array([1, 2, 3, 4, 1, 2, 3, 4, 1, 2], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings(ind, 10, 7)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 6)
        self.assertEqual(cvnseg, 4)

        # A = 10, K = 4
        ind = np.array([1, 2, 3, 4, 1, 2, 3, 4, 1, 2], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings(ind, 10, 10)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 6)
        self.assertEqual(cvnseg, 4)


    def test_get_settings_random(self):
        """
        Test that "get_settings" return correct settings for random splits.
        """

        # A = 5, K = 3
        ind = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'rand', 'nseg': 3}, 10, 5)
        np.testing.assert_array_equal(np.sort(cvind), np.sort(ind))
        self.assertRaises(AssertionError, np.testing.assert_array_equal, cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 10, K = 4
        ind = np.array([1, 2, 3, 4, 1, 2, 3, 4, 1, 2], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'rand', 'nseg': 4}, 10, 10)
        np.testing.assert_array_equal(np.sort(cvind), np.sort(ind))
        self.assertRaises(AssertionError, np.testing.assert_array_equal, cvind, ind)
        self.assertEqual(cvncomp, 6)
        self.assertEqual(cvnseg, 4)

        # A = 10, K = 10
        ind = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'rand', 'nseg': 10}, 10, 10)
        np.testing.assert_array_equal(np.sort(cvind), np.sort(ind))
        self.assertRaises(AssertionError, np.testing.assert_array_equal, cvind, ind)
        self.assertEqual(cvncomp, 8)
        self.assertEqual(cvnseg, 10)


    def test_get_settings_ven(self):
        """
        Test that "get_settings" return correct settings for venetian blinds splits.
        """

        resp = np.array([9.1, -8.0, 1.0, 2.0, 3.1, -11.5, 1.2, 4.5, 9.1, -2.2])

        # information for help
        #  resp_sorted = [-11.5, -8.0, -2.2, 1.0, 1.2, 2.0, 3.1, 4.5, 9.1, 9.1]
        #  rind1 = [5, 1, 9, 2, 6, 3, 4, 7, 0, 8]
        #  rind2 = [8, 1, 3, 5, 6, 0, 4, 7, 9, 2]
        #  ind_sorted = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1]
        #  ind_unsorted = [3, 2, 1, 3, 1, 1, 2, 2, 1, 3]
        #  seg_sorted = (1: -11.5, 1.0, 3.1, 9.1; 2: -8.0, 1.2, 4.5; 3: -2.2, 2.2, 9.1)
        #  seg_unsorted = (1: 9.1, 1.0, 3.1, 11.5; 2: -8.0, 1.2, 4.5; 3: -2.2, 2.2, 9.1)


        # A = 5, K = 3, sorted resp
        ind = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 3}, 10, 5, np.sort(resp))
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 5, K = 3, unsorted resp
        ind = np.array([3, 2, 1, 3, 1, 1, 2, 2, 1, 3], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 3}, 10, 5, resp)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 5, K = 3, no resp
        ind = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 3}, 10, 5)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 5)
        self.assertEqual(cvnseg, 3)

        # A = 10, K = 10, sorted resp
        ind = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 10}, 10, 10, np.sort(resp))
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 8)
        self.assertEqual(cvnseg, 10)

        # A = 10, K = 10, unsorted resp
        ind = np.array([9, 2, 4, 6, 7, 1, 5, 8, 10, 3], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 10}, 10, 10, resp)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 8)
        self.assertEqual(cvnseg, 10)

        # A = 10, K = 10, no resp
        ind = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype = int)
        (cvind, cvncomp, cvnseg) = get_cvsettings({'type': 'ven', 'nseg': 10}, 10, 10)
        np.testing.assert_array_equal(cvind, ind)
        self.assertEqual(cvncomp, 8)
        self.assertEqual(cvnseg, 10)


if __name__ == '__main__':
    unittest.main()