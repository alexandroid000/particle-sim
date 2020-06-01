import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from utilities import *
from configuration import *
from copy import copy
from time import sleep


policy = []
with open("policy9.txt", 'r') as p:
    line = p.readline().strip().strip('()').split(", ")
    policy = [int(p) for p in line]

N = len(policy)

counts = {}

for p in policy:
    if p in counts:
        counts[p] += 1
    else:
        counts[p] = 1

print(len(counts), "actions used")
fractions = [(k, v/N) for (k,v) in counts.items()]
print(fractions)

policy25 = ['CCW', 'X', 'X', 'CW']
wire_verts = np.array(mk_regpoly(4, 0.4*L, offset=np.pi/4))
wires = [Wire(v, o) for v, o in zip(wire_verts, policy25)]

# Grid of x, y points
nx, ny = 10, 10 
x = np.linspace(-L, L, nx)
y = np.linspace(-L, L, ny)
X, Y = np.meshgrid(x, y)


Bx, By = np.zeros((ny, nx)), np.zeros((ny, nx))
i = 0
for wire in wires:
    Xi, Yi = copy(X), copy(Y)
    bx, by = wire.force_at(Xi, Yi)
    Bx += bx
    By += by
    i += 1


fig = plt.figure()
ax = fig.add_subplot(111)

plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')

ax.streamplot(x, y, Bx, By, linewidth=1, density=1, arrowstyle='->', arrowsize=1.5)

#plt.quiver(X, Y, Bx, By)
for w in wires:
    pos = w.xy
    c = Circle(pos, 0.1, facecolor='r')
    ax.add_patch(c)


oct_verts = np.array(mk_regpoly(8, 0.8*L, offset=np.pi/8.))
p = Polygon(oct_verts, ec='k', lw=2, fc='none')
ax.add_patch(p)

for r in rs_as_obs[1:]:
    p = Polygon(r, ec='k', lw=2, linestyle=':', fc='none')
    ax.add_patch(p)

R = 0.85*L

ax.set_xlim(-R,R)
ax.set_ylim(-R,R)
ax.set_aspect('equal')
plt.savefig("policy25field.eps", bbox='tight')

times = [4.27, 21.72, 119.4, 644.17, 3413.6, 18626.6, 98147.5, 522008.7]

time_dat = zip(range(2,10), times)

plt.clf()
ax = fig.add_subplot(111)
plt.plot(range(2,10), times, marker='o', linewidth=2, markersize=12)
plt.xlabel("Number of Agents")
plt.ylabel("Time (s)")
plt.title("Time to Converge MDP Policy with Value Iteration")
plt.savefig("MDPtiming.eps", bbox='tight')
