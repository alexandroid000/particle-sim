#! /usr/bin/env python

import unittest
from backend import *
from utilities import *
from random import randint, random
import numpy as np


class TestExamples(unittest.TestCase):

    # add global stuff here
    def setUp(self):
        return

    def test_scatter(self):
        p1 = [0.,0.]
        p2 = [1.,0.]
        m1, m2 = 1., 1.
        v1 = [1.,1.]
        v2 = [-1., 1.]

        particle1 = Particle(position= p1,
                             velocity= v1,
                             radius = None,
                             species= 'A-free',
                             mass = m1)

        particle2 = Particle(position= p2,
                             velocity= v2,
                             radius = None,
                             species= 'A-free',
                             mass = m2)

        twoCollide(particle1, particle2)

        np.testing.assert_almost_equal(particle1.velocity,
                                       [-1., 1.],
                                       decimal=7, verbose=True)
        np.testing.assert_almost_equal(particle2.velocity,
                                       [1., 1.],
                                       decimal=7, verbose=True)


    def test_magnet(self):
        o = np.array([0.,0.])
        pt1 = np.array([random(),random()])
        pt2 = o + 2*pt1

        wire1 = Wire(o, "CCW")
        wire2 = Wire(o, "CW")
        pt1_test_ccw = np.array(wire1.force_at(*pt1))
        pt2_test_ccw = np.array(wire1.force_at(*pt2))
        pt1_test_cw  = np.array(wire2.force_at(*pt1))
        pt2_test_cw  = np.array(wire2.force_at(*pt2))
        np.testing.assert_almost_equal(pt1_test_ccw,
                                       rotate_vector(pt1_test_cw, np.pi),
                                       decimal=7, verbose=True)
        np.testing.assert_almost_equal(pt1_test_ccw/2.,
                                       pt2_test_ccw,
                                       decimal=7, verbose=True)

    def test_bb(self):
        oct = Simple_Polygon("oct",np.array(mk_regpoly(8, 5.)))
        min_x, max_x, min_y, max_y, bb = mk_bounding_box(oct)
        np.testing.assert_almost_equal(min_x, -5., decimal=7, verbose=True)
        np.testing.assert_almost_equal(max_x, 5., decimal=7, verbose=True)
        np.testing.assert_almost_equal(min_y, -5., decimal=7, verbose=True)
        np.testing.assert_almost_equal(max_y, 5., decimal=7, verbose=True)

    def test_encode_decode(self):
        test1 = ["CCW", "X", "CW", "CW"]
        test2 = ["X", "X", "X", "X"]
        test3 = ["CCW", "CCW", "CCW", "CCW"]

        self.assertSequenceEqual(decode_policy(encode_policy(test1)), test1)
        self.assertSequenceEqual(decode_policy(encode_policy(test2)), test2)
        self.assertSequenceEqual(decode_policy(encode_policy(test3)), test3)

