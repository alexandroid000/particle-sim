from utilities import *
from environments import *
import numpy as np

# System Configuration
# --------------------

# define simulation parameters here

L = 150.0 # environment characteristic length
N = 9 # number of particles
T = 300 # number of time steps/stages in simulator, unitless 
R = 1 # colliding length scale ## Any units? 1/50 of animation scale?
border_region = R
allow_attachment = False # particles don't interact with each other #
xMin = -2*L
xMax = 2*L
yMin = -2*L
yMax = 2*L

## play around with numbers below 

# type A particles:
    # faster
    # smaller rotational drift
    # escape from walls more quickly
A_properties = {'vel':25.0, 'wall_prob': 0.05, 'beta': 0.2, 'mass':1}

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

env = simple_nonconv_p()

