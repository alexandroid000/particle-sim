#! /usr/bin/env python

import unittest
from backend import *
from utilities import *
from environments import *
from random import randint, random
import numpy as np


class TestExamples(unittest.TestCase):

    # add global stuff here
    def setUp(self):
        return

    def test_grid(self):
        system = System()
        L = 1.0
        T = 10
        N = 4
        R = 0.25
        BR = 0.005
        K = 0.5
        ATTACH = False
        data = {"pos":[[]]*T, "env":[[]]*T, "counts":[[]]*(T-1), "wires":[[]]*T}
        env = square(L)

        start_pts = [np.array([0.25, 0.25]),
                     np.array([0.25, 0.75]),
                     np.array([0.75, 0.25]),
                     np.array([0.75, 0.75])]

        for i in range(N):
            vel = normalize(np.array([random()-0.5, random()-0.5]))
            system.particles.append(Particle(position=start_pts[i],
                                            velocity=list(vel),
                                            radius = R,
                                            species= 'A-free',
                                            mass = 100.0))

        simulation = ParticleSim(system, data, env,
                                 br = BR, k = K, sticky=ATTACH,
                                 r = R)

        Ncells = 0
        for f in list(simulation.cells.items()):
            print(f)
        for x,ys in simulation.cells.items():
            for y,c in ys.items():
                Ncells += 1

        self.assertEqual(Ncells, 4)

        np.testing.assert_almost_equal(simulation.cells[0][0].position, start_pts[0], decimal=7, verbose=True)
        np.testing.assert_almost_equal(simulation.cells[0][1].position, start_pts[1], decimal=7, verbose=True)
        np.testing.assert_almost_equal(simulation.cells[1][0].position, start_pts[2], decimal=7, verbose=True)
        np.testing.assert_almost_equal(simulation.cells[1][1].position, start_pts[3], decimal=7, verbose=True)


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

