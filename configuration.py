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
BR: float
    characteristic "collision" length scale with environment boundaries
K: float
    magnitude of repulsive force between particles
ATTACH: bool
    true if we want particles to stick together when they collide
ANIMATE: bool
    true if we want to produce a mp4 of the simulation
FRAC_BALLISTIC: float
    should be between 0 and 1. fraction of particles of type B, ballistic
"""

L = 1.0
R = 0.05*L
BR = 0.005
K = 0.5
ATTACH = False
ANIMATE = True

# following are set via command line or defaults now

#N = 30
#T = 100
#FRAC_BALLISTIC = 0.2

"""env defines the environment

needs to be an instance of class
Simple_Polygon from bounce-viz/src/simple_polygon.py

see environments.py and bounce-viz/src/maps.py for more examples"""
#env = Simple_Polygon("smallpoly1",smallpoly1[0])
env = octagon(L)

"""bounds of animation window"""
xs = [x for (x,y) in env.complete_vertex_list]
ys = [y for (x,y) in env.complete_vertex_list]
XMIN = np.amin(xs)
XMAX = np.amax(xs)
YMIN = np.amin(ys)
YMAX = np.amax(ys)



# Particle Configuration
# ----------------------

## play around with numbers below 

# type A particles:
    # active brownian particles
A_properties = {'vel':1., 'wall_prob': 0.05, 'beta': 0.2}

# type B particles:
    # bouncing robots 
B_properties = {'vel':200., 'wall_prob': 0.05, 'beta': 0.0}

properties = { 'A-free':A_properties
             ,'B-free':B_properties
             ,'A-wall':A_properties
             ,'B-wall':B_properties
             }
