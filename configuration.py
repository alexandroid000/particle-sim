from utilities import *
from environments import *
import numpy as np

# System Configuration
# --------------------

# define simulation parameters here

L = 3.0 # environment characteristic length
N = 9 # number of particles
T = 300 # number of time steps/stages in simulator, unitless 
R = 0.02 # colliding length scale; effective radius of particle
ATTACH = True # whether particles stick to each other
XMIN = -2*L
XMAX = 2*L
YMIN = -2*L
YMAX = 2*L
ANIMATE = True

# Environment Configuration
# --------------------

env = square(L)

# Particle Configuration
# ----------------------

# type A particles:
    # faster
    # smaller rotational drift
    # escape from walls more quickly
A_properties = {'vel':1.0, 'wall_prob': 0.05, 'beta': 0.2}

# type B particles:
    # slower
    # more rotational drift
    # get stuck on walls more
B_properties = {'vel':0.3, 'wall_prob': 0.2, 'beta': 0.5}

properties = { 'A-free':A_properties
             ,'B-free':B_properties
             ,'A-wall':A_properties
             ,'B-wall':B_properties
             }

