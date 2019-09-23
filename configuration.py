from utilities import *
from environments import *
import numpy as np

# System Configuration
# --------------------

"""
Parameters
----------

L: float
    a characteristic environment length, used to autogenerate some environments
N: int
    the number of particles at the beginning of the simulation
T: int
    number of time steps in simulator
R: float
    characteristic "collision" length scale: width of bounding box for particles
ATTACH: bool
    true if we want particles to stick together when they collide
ANIMATE: bool
    true if we want to produce a mp4 of the simulation
"""

L = 3.0
N = 200
T = 300
R = 0.02
ATTACH = False
ANIMATE = True

"""env defines the environment

needs to be an instance of class
Simple_Polygon from bounce-viz/src/simple_polygon.py

see environments.py and bounce-viz/src/maps.py for more examples"""
env = Simple_Polygon("poly1",poly1[0])

"""bounds of animation window"""
xs = [x for (x,y) in env.complete_vertex_list]
ys = [y for (x,y) in env.complete_vertex_list]
XMIN = np.amin(xs)
XMAX = np.amax(xs)
YMIN = np.amin(ys)
YMAX = np.amax(ys)



# Particle Configuration
# ----------------------

# type A particles:
    # faster
    # smaller rotational drift (beta)
    # escape from walls more quickly (wall_prob, probability of staying stuck on wall)
A_properties = {'vel':20.0, 'wall_prob': 0.05, 'beta': 0.2}

# type B particles:
    # slower
    # more rotational drift
    # get stuck on walls more
B_properties = {'vel':7.0, 'wall_prob': 0.2, 'beta': 0.5}

properties = { 'A-free':A_properties
             ,'B-free':B_properties
             ,'A-wall':A_properties
             ,'B-wall':B_properties
             }

