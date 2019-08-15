from random import random
import sys

sys.path.append('../')
sys.path.append('../bounce-viz/src/')

from backend import *
from utilities import *


# type of v1 is two-vector
def twoCollide(m1, v1, m2, v2):
    v1x, v1y = v1
    v1prime = v1
    v2prime = v2
    return v1prime, v2prime

v1 = normalize(np.array([1.0, 0.0]))
v2 = normalize(np.array([-1.0, 0.0]))

m1 = 1.0
m2 = 2.0

if __name__ == '__main__':
    print("computed v1prime and v2prime are:")
    v1prime, v2prime = twoCollide(m1, v1, m2, v2)
    print(v1prime, v2prime)
