#! /usr/bin/env python

import unittest
from backend import *
from utilities import *
from maps import *
from environments import *
from random import randint, random
import numpy as np


class TestExamples(unittest.TestCase):

    # add global stuff here
    def setUp(self):
        return

    def test_discretize(self):

        rect = np.array([(0., 0.), (12., 0.), (12., 8.), (0., 8.)])

        env = Simple_Polygon("r", rect)
        d = discretizedEnvironment(env, N=12)
        N_x = len(d.cells)
        N_y = len(d.cells[0])
        self.assertEqual(N_x, 12)
        self.assertEqual(N_y, 8)

        particle_pos = []

        p1 = (0.5,0.5)
        p2 = (-1.,-1.)
        p3 = (13., 13.)
        p4 = (5.5, 5.5)

        self.assertEqual(d.quadrant(*p1), (0,0))
        self.assertRaises(ValueError, d.quadrant, *p2)
        self.assertRaises(ValueError, d.quadrant, *p3)
        self.assertEqual(d.quadrant(*p4), (5,5))


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


    def test_repulse(self):
        p1 = [0.,0.]
        p2 = [1.,0.]
        m1, m2 = 1., 1.
        v1 = [1.,0.]
        v2 = [-1., 0.]

        particle1 = Particle(position= p1,
                             velocity= v1,
                             radius = 1.0,
                             species= 'A-free',
                             mass = m1)

        particle2 = Particle(position= p2,
                             velocity= v2,
                             radius = 1.0,
                             species= 'A-free',
                             mass = m2)

        softRepulse(particle1, particle2, 1)

        np.testing.assert_almost_equal(particle1.velocity,
                                       [0., 0.],
                                       decimal=7, verbose=True)
        np.testing.assert_almost_equal(particle2.velocity,
                                       [0., 0.],
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

