from utilities import *
from environments import *
import numpy as np

# System Configuration
# --------------------

# define simulation parameters here

L = 3.0 # environment characteristic length
N = 9 # number of particles
T = 300 # number of time steps/stages in simulator, unitless 
R = 0.02 # colliding length scale 
border_region = R
allow_attachment = False # particles don't interact with each other

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

# Environment Configuration
# --------------------

env = square(L)

