from random import random
import sys
import numpy as np
sys.path.append('../')
sys.path.append('../bounce-viz/src/')

from backend import *
from utilities import *


# type of v1 is two-vector: 1-D 
## atan2 use this for 2-d
###graph grammer 
#mi : mass vi: velocity pi: position
v1 = np.array([1.0, -3.0]) #initial speed # WHY normalize? deleted normalization and it worked! 
v2 = np.array([4.0, 8.0])  # 4 8 
p1 = np.array([1,2])
p2 = np.array([1,0])
m1 = 1.0
m2 = 2.0
def twoCollide(m1, v1, m2, v2,p1,p2):
    v1x, v1y = v1
    v2x, v2y = v2
    x1x, x1y = p1 #p is vector of particle 1's position
    x2x, x2y = p2
    v1prime = v1 - (2 * m2) / (m1 + m2) * (np.dot(v1-v2,p1-p2)) / ((x1x - x2x)**2 + (x1y - x2y)**2)  * (p1 - p2)
    v2prime = v2 - (2 * m1) / (m1 + m2) * (np.dot(v2-v1,p2-p1)) / ((x2x - x1x)**2 + (x2y - x1y)**2) * (p2 - p1)
    return v1prime, v2prime
## Test with vector calculation. test for 2D


if __name__ == '__main__':
    print("computed v1prime and v2prime are:")
    v1prime, v2prime = twoCollide(m1, v1, m2, v2, p1, p2)
    print(v1prime, v2prime)


#def ncollide(m1, v1, M2, V2 ): # n nuber of particles # stick or bounce off # more than 2 particles
    #v1x 