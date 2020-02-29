#! /usr/bin/env python

import csv
from random import random

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
import matplotlib.animation as animation

from backend import *
from configuration import *
from analyze_density import get_limit_cycle_poly
from utilities import normalize, uniform_sample_along_circle, uniform_sample_from_poly

# Animation display parameters
# ----------------------------

color_map = { 'A-wall': (0, 0, 1)
            , 'A-free': (0, 0, 1)
            , 'B-wall': (0, 1, 0)
            , 'B-free': (0, 1, 0)
            }
size_map =  { 'A-wall': 10
            , 'A-free': 10
            , 'B-wall': 50
            , 'B-free': 50
            }

# clean data and put in form that makes animate functions more readable
  
class Data:

    def __init__(self, db, start=0):
        self.xy = db["pos"]
        self.env = db["env"]
        self.num = start

    def __iter__(self):
        return self

    def clean_system(self, dat):
        types = [t for (t,xy) in dat]
        xys = [xy for (t,xy) in dat]
        return types, np.array(xys)

    def __next__(self):
        dat = self.xy[self.num]
        polys = [np.array(poly) for poly in self.env[self.num]]
        self.num += 1
        return self.clean_system(dat), polys


def write_data(database, simname):
        
    # write position data to file
    # only log type A particles
    with open(simname+'_typeA.xyz','w') as th:
        for i in range(T):
            xys = database["pos"][i]
            for (t, [x,y]) in xys:
                if t[0] == 'A':
                    th.write(str(x)+" "+str(y)+" ")
            th.write("\n")

    print("wrote data to",simname+"_typeA.xyz")

def init():
    """initialize animation"""
    global scat, patches
    patches = []

    for poly in initenv:
        p = Polygon(poly, ec='k', lw=2, fc='none')
        patches.append(ax.add_patch(p))

    cycle = get_limit_cycle_poly(8, 1., 1, 0.2)
    c = Polygon(cycle, ec='r', lw=1, fc='none', ls='--')
    cpatch = [ax.add_patch(c)]

    return patches+cpatch+[scat]

def animate(i):
    """perform animation step"""
    global scat, patches
    xy, polys = next(d)

    for j in range(1,len(patches)):
        patches[j].set_xy(polys[j])
        patches[j].set_edgecolor('k')

    # update pieces of the animation
    scat.set_facecolors([color_map[t] for t in xy[0]])
    scat.set_sizes([size_map[t] for t in xy[0]])
    scat.set_offsets(xy[1])
    return patches+[scat]

def mkAnimation():
    ani = animation.FuncAnimation(fig, animate, frames=T-2, interval=10,
                              blit=True, init_func=init)
    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    ani.save(simname+'.mp4', fps=15, extra_args=['-vcodec', 'libx264'])


if __name__ == '__main__':

    # TODO: use a typed arguments library here with default values etc
    args = sys.argv[1:]
    if len(args) == 0:
        FRAC_BALLISTIC = 0.1
        N = 300
        T = 400
    elif len(args) == 1:
        FRAC_BALLISTIC = float(args[0])*0.1
        N = 300
        T = 10
    elif len(args) == 2:
        FRAC_BALLISTIC = float(args[0])
        N = float(args[1])
        T = 400
    else:
        FRAC_BALLISTIC = float(args[0])
        N = float(args[1])
        T = float(args[2])


    # initialize simulation
    system = System()
    data = {"pos":[[]]*T, "env":[[]]*T}
    simulation = ParticleSim(system, data, env,
                             br = BR, k = K, sticky=ATTACH,
                             r = R)
    simname = env.name+"_N"+str(N)+"_T"+str(T)+"_F"+str(FRAC_BALLISTIC)

    # create N particles at random locations in the polygon
    start_pts = uniform_sample_from_poly(env, N)
    #start_pts = uniform_sample_along_circle(env, N, 2.0)
    for i in range(N):
        vel = normalize(np.array([random()-0.5, random()-0.5]))
        if random() < FRAC_BALLISTIC:
            system.particle.append(Particle(position=start_pts[i],
                                            velocity=list(vel),
                                            radius = R,
                                            species= 'B-free',
                                            mass = 100.0))

        else:
            system.particle.append(Particle(position=start_pts[i],
                                            velocity=list(vel),
                                            radius = R,
                                            species= 'A-free',
                                            mass = 1.0))
    # run simulation for T steps
    simulation.run(T-1)
    print("Finished Simulation, Writing Data...")

    write_data(simulation.db, simname)

    if ANIMATE:
        # make iterator to make animation easier
        d = Data(simulation.db)
        initxy, initenv = copy(next(d))

        colors = [color_map[t] for t in initxy[0]]
        sizes = [size_map[t] for t in initxy[0]]

        fig = plt.figure()
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                             xlim=(XMIN,XMAX) , ylim=(YMIN,YMAX)) #scale for animation window
        
        scat = ax.scatter(initxy[1][:,0]
                        , initxy[1][:,1]
                        , facecolors = colors
                        , s = sizes
                        )

        print("writing video to",simname+".mp4")
        mkAnimation()


