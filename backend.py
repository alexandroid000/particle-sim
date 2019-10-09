import numpy as np
from copy import copy, deepcopy
from random import random


# using bounce-viz as a submodule for geometric utilities
import sys
sys.path.insert(0, "./bounce-viz/src/")
from helper.shoot_ray_helper import IsInPoly, ClosestPtAlongRay # pylint: disable=unused-import
from helper.geometry_helper import AngleBetween # pylint: disable=unused-import
from utilities import twoCollide
from configuration import * # pylint: disable=unused-import

###IDEA: Mass--> different types of collision ( in run function), different scatter
# Simulation Backend
# ------------------

class Particle():

    def __init__(self, position, velocity, radius = None, species = None, mass = None): # constructoe of  a class
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.radius = radius
        self.species = species
        self.mass = mass
### Q. why position and velocity are not set equal to "None"? Design choice? 
## self is a convention to write 
## Mass should be defined in Species of Particles 
class System(): #list of particles

    def __init__(self):
        self.particle = []

class ParticlePhysics(object):

    def __init__(self, system, env, delta=0.05,
                       br = 0.01, k = 1.0, sticky = True,
                       wires = []):
        self.system = system
        self.env = env
        self.delta = delta
        self.vs = self.env.vertex_list_per_poly
        self.n = len(self.vs[0]) # outer boundary
        self.br = br
        self.K = k
        self.sticky = sticky
        self.wires = wires

    def neighbors(self, particle): #distinguishes if particles are neighboring each other 
        [x,y] = particle.position
        neighbors = []
        bounding_box = Simple_Polygon("bb",np.array([[x-R,y-R]
                                                   , [x+R,y-R]
                                                   , [x+R,y+R]
                                                   , [x-R,y+R]]))
        for p in self.system.particle:
            if IsInPoly(p.position, bounding_box) and (p is not particle): # does not count current particle as its own neighbor 
                neighbors.append(p)                                        # confusing variable names
        return neighbors
    # read about scoping 
   
    # TODO: replace this with polygon offset calculator to make more robust for
    # nonconvex polygons
    # look into pyclipper library
    def project_to_border_region(self, pt):
        d, [ex,ey] = dist_dir_closest_edge(pt, self.env)
        normal_dir = normalize(np.array([-ey, ex])) # pointing into polygon
        return pt + (self.br - d) * normal_dir

    # move obstacle #c in the direction of dir
    # currently only internal obstacles can move, boundary is fixed
    def move_obstacle(self, c, dir):
        old_poly = [[v for (i,v) in c] for c in deepcopy(self.env.vertex_list_per_poly)]
        new_poly = [(v + dir) for v in old_poly[c]]
        new_obstacles = old_poly[:c]+[new_poly]+old_poly[(c+1):]
        print("moved obstacle",c,"along vector",dir)
        print(len(new_obstacles),"new obstacles:",new_obstacles)
        new_env = Simple_Polygon(self.env.name, new_obstacles[0], new_obstacles[1:])
        self.env = new_env

    # only called if a moveable obstacle (component c != 0)
    def obstacle_interaction(self, particle, edge_dir):
        d, c, j = closest_edge(particle.position, self.env)
        if c != 0:
            dr = self.next_dr(particle) #picks new direction
            [ex,ey] = edge_dir
            proj_dr_onto_edge = (dr.dot(edge_dir))/np.linalg.norm(edge_dir)
            particle.position += proj_dr_onto_edge

            push_dir = normalize(np.array([ey, -ex])) # pointing into obstacle
            proj_dr_onto_normal = (dr.dot(push_dir))/np.linalg.norm(push_dir)
            self.move_obstacle(c, proj_dr_onto_normal)


    # how the particles should leave the static obstacles
    def scatter(self, particle, edge_dir):
        particle.species = particle.species[0]+'-free'
        [ex,ey] = edge_dir
        normal = normalize(np.array([-ey, ex])) # pointing into polygon

        # rotate particle's velocity uniformly out from wall
        th_out = np.pi*random() - np.pi/2
        particle.velocity = normalize(rotate_vector(normal, th_out))

    # TODO: order of collisions is arbitrary if more than one neighbor in
    # bounding box. Is there a more principled way to do this?
    def group_collision(self, particle, ns):
        pairs = [(particle, n) for n in ns]
        for (particle,n) in pairs:
            #twoCollide(particle, n) # uncomment for elastic collision
            softRepulse(particle, n, self.K)

    def take_step_boundary(self, particle):
        d, edge_dir = dist_dir_closest_edge(particle.position, self.env)

        # escape from wall
        if random() > properties[particle.species]['wall_prob']:
            self.scatter(particle, edge_dir)

        # stay on wall, impart force
        else:
            self.obstacle_interaction(particle, edge_dir)

    def next_dr(self, particle): #Next direction 

        # Brownian motion, random step on unit circle
        [xi_x, xi_y] = normalize([random()-0.5, random()-0.5])

        # compute direction of next step
        v = properties[particle.species]['vel']
        xdot = v*particle.velocity[0] + xi_x
        ydot = v*particle.velocity[1] + xi_y

        # random update to velocity heading
        # Gaussian, mean at current heading, standard deviation 1 radian
        theta = np.arctan2(particle.velocity[1], particle.velocity[0])
        xi_theta = np.random.vonmises(mu=theta, kappa=1.)

        # update velocity for next step
        beta = properties[particle.species]['beta']
        if self.wires != []: 
            particle.velocity = force_from_wires(self.wires, particle.position)
        particle.velocity[0] += beta*np.cos(xi_theta)
        particle.velocity[1] += beta*np.sin(xi_theta)
        particle.velocity = normalize(particle.velocity)

        # take step, scaled by delta
        dr = self.delta*np.array([xdot, ydot])
        return dr

    def take_step(self, particle):
        dr = self.next_dr(particle)
        if IsInPoly(particle.position + dr, self.env):
            particle.position += dr
        else: #prevents particles going outside of polygon 
            particle.position = self.project_to_border_region(particle.position)
            particle.species = particle.species[0]+'-wall'

class ParticleSim(ParticlePhysics):

    def __init__(self, system, database, env, delta=0.02,
                       br = 0.01, k = 1.0, sticky = True, wires = [],
                       regions = [], policy = []):

        ParticlePhysics.__init__(self, system, env, delta, br, sticky, wires)

        self.system = system #list of particle
        self.db = database
        self.delta = delta
        self.env = env
        self.vs = self.env.vertex_list_per_poly
        self.n = len(self.vs[0]) # outer boundary
        self.br = br
        self.K = k
        self.sticky = sticky
        self.wires = wires
        self.regions = regions
        self.policy = policy


    def particle_collide(self, p): 
        ns = self.neighbors(p) # checks bounding box for neighbors
        if self.sticky and ns != []:
            # TODO: add probability of sticking
            # eventually, dependent on shape/angle of incidence
            mode = p.species[2:]
            p.species = 'B-'+ mode
            # each particle is an object, when colliding all but one get removed
            # TODO: update mass of the mega-particle
            for n in ns:
                self.system.particle.remove(n)
        else:
            if ns != []:
                self.group_collision(p, ns)


    def run(self, steps):
        print("Running Simulation for",steps,"steps")

        # initialize region counts
        region_counts = [0]*len(self.regions)
        states = [0]*len(self.system.particle)
        for j,p in enumerate(self.system.particle):
            for i,r in enumerate(self.regions):
                if IsInPoly(p.position, r):
                    region_counts[i] += 1
                    states[j] = i

        # run sim for T steps
        ##board pic 
        for i in range(steps):

            if i % 10 == 0:
                print("Step",i)

            # log regions; only works at beginning of loop for some reason
            joint_state = encodeJointState(states)
            if self.policy != []:
                new_orientations = decode_policy(self.policy[joint_state])
                self.log_data(i, region_counts, new_orientations)
                self.wires = [Wire(v, o) for v, o in zip(wire_verts, new_orientations)]
            else:
                self.log_data(i, region_counts)

            region_counts = [0]*len(self.regions)
            for j,p in enumerate(self.system.particle):
                for i,r in enumerate(self.regions):
                    if IsInPoly(p.position, r):
                        region_counts[i] += 1
                        states[j] = i

                # detect particle-particle collisions
                self.particle_collide(p)
                # boundary mode
                if p.species[2:] == "wall":
                    self.take_step_boundary(p)
                # free space mode
                else:
                    self.take_step(p)


    def log_data(self, step, r_counts, wires = []):
        xys = [(copy(p.species), copy(p.position)) for p in self.system.particle]
        envs = [[v for (i,v) in c] for c in deepcopy(self.env.vertex_list_per_poly)]
        rs = copy(r_counts)
        wires = deepcopy(wires)
        self.db["pos"][step] = xys
        self.db["env"][step] = envs
        self.db["counts"][step] = rs
        self.db["wires"][step] = wires

