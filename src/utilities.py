import numpy as np
import sys
sys.path.insert(0, "./bounce-viz/src/")
from simple_polygon import Simple_Polygon
from random import uniform
from helper.shoot_ray_helper import IsInPoly, ClosestPtAlongRay

EPSILON = 0.000000001

# Geometric Operations
# --------------------

# http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
def closest_edge(pt, poly):
    [x0,y0] = pt
    vs = poly.vertex_list_per_poly
    n = poly.size
    components = len(vs)
    min_d = 100000000000
    closest_component = -1
    closest_edge = -1
    # find closest edge over external boundary and holes
    for (i, component) in enumerate(vs):
        m = len(component)
        for j in range(m):
            [x1, y1], [x2,y2] = component[j][1], component[(j+1) % m][1]
            d = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / np.sqrt((x2-x1)**2 + (y2-y1)**2)
            if d < min_d:
                min_d = d
                closest_component = i
                closest_edge = j
    return min_d, closest_component, closest_edge

def dist_dir_closest_edge(pt, poly):
    d, c, j = closest_edge(pt, poly)
    vs = poly.vertex_list_per_poly
    csize = len(vs[c])
    edge_vect = vs[c][(j + 1) % csize][1] - vs[c][j][1]
    return d, edge_vect

def normalize(vector):
    norm = np.linalg.norm(vector)
    if norm > EPSILON:
        return vector/norm
    else:
        return 0.0*vector

def rotate_vector(v, theta):
    vx, vy = v[0], v[1]
    return np.array( [np.cos(theta)*vx - np.sin(theta)*vy,
                     np.sin(theta)*vx + np.cos(theta)*vy])

def midpoint(pt1, pt2):
    return (pt1+pt2)/2


# Collision Utilities
# -------------------


# elastic scattering
def twoCollide(particle1, particle2):
    v1, v2 = particle1.velocity, particle2.velocity
    v1x, v1y = v1
    v2x, v2y = v2
    m1, m2 = particle1.mass, particle2.mass
    p1, p2 = particle1.position, particle2.position
    x1x, x1y = p1
    x2x, x2y = p2
    v1prime = v1 - (2 * m2) / (m1 + m2) * (np.dot(v1-v2,p1-p2)) / ((x1x - x2x)**2 + (x1y - x2y)**2)  * (p1 - p2)
    v2prime = v2 - (2 * m1) / (m1 + m2) * (np.dot(v2-v1,p2-p1)) / ((x2x - x1x)**2 + (x2y - x1y)**2) * (p2 - p1)
    particle1.velocity = v1prime
    particle2.velocity = v2prime

def softRepulse(particle1, particle2, K):
    v1, v2 = particle1.velocity, particle2.velocity
    v1x, v1y = v1
    v2x, v2y = v2
    m1, m2 = particle1.mass, particle2.mass
    p1, p2 = np.array(particle1.position), np.array(particle2.position)
    r1, r2 = particle1.radius, particle2.radius
    # vector from p1 to p2
    r12 = normalize(p2-p1)
    # vector from p2 to p1
    r21 = normalize(p1-p2)

    # magnitude of repulsive force proportional to distance
    rlen = np.linalg.norm(p1-p2)
    if rlen < r1+r2:
        f = K*(r1 + r2 - rlen)
    else:
        f = 0

    # force on particle 1
    # TODO: fix to incorporate direction of heading
    f1 = f*r21
    # force on particle 2
    f2 = f*r12

    particle1.velocity += f1/m1
    particle2.velocity += f2/m2

class discretizedEnvironment():
    def __init__(self, env, D):
        self.env = env
        self.D = D # number of cells along longest axis
        xs = [x for (x,y) in env.complete_vertex_list]
        ys = [y for (x,y) in env.complete_vertex_list]
        self.XMIN = np.amin(xs)
        self.XMAX = np.amax(xs)
        self.YMIN = np.amin(ys)
        self.YMAX = np.amax(ys)

        L_actual = max(abs(self.XMAX-self.XMIN), abs(self.YMAX-self.YMIN))
        self.R = L_actual / self.D
        xnum = int(np.ceil(abs(self.XMAX-self.XMIN)/self.R))
        ynum = int(np.ceil(abs(self.YMAX-self.YMIN)/self.R))
        self.cells = {i:{j:[] for j in range(ynum)} for i in range(xnum)}
        self.L = len(self.cells) # number of rows
        self.M = len(self.cells[0]) # number of columns

    def insert_particle(self, id, x, y):
        (i,j) = self.quadrant(x,y)
        self.cells[i][j].append(id)

    def remove_particle(self, id, i, j):
        self.cells[i][j].remove(id)

    def quadrant(self, x,y):
        if x > self.XMAX or x < self.XMIN or y > self.YMAX or y < self.YMIN:
            raise ValueError("Particle has escaped the polygon...")
        r_num_x = abs(x-self.XMIN) // self.R
        r_num_y = abs(y-self.YMIN) // self.R
        return (int(r_num_x), int(r_num_y))

    def update(self, particle, old_i, old_j):
        id = particle.id
        [x,y] = particle.position
        i,j = self.quadrant(x,y)
        if (i != old_i) or (j != old_j):
            self.remove_particle(id, old_i, old_j)
            self.insert_particle(id, x, y)



# Environment Utilities
# -----------------

def mk_spiky_circle(n, r):
    d = 2*np.pi/n
    theta = 0
    pts = []
    for i in range(n):
        pt1 = [r*np.cos(theta), r*np.sin(theta)]
        r2 = 1.5*r
        pt2 = [r2*np.cos(theta), r2*np.sin(theta)]
        theta += d
        pts.extend([pt1, pt2])
    return pts

def mk_spiky_obstacle(n, r):
    d = 2*np.pi/n
    theta = 0.0
    pts = []
    for i in range(n):
        pt1 = [r*np.cos(theta), r*np.sin(theta)]
        r2 = 0.6*r
        th2 = theta + 0.3*d
        pt2 = [r2*np.cos(th2), r2*np.sin(th2)]
        theta += d
        pts.extend([pt1, pt2])
    return pts[::-1]

def mk_regpoly(n, r, offset=0.0):
    d = 2*np.pi/n
    theta = offset
    pts = []
    for i in range(n):
        pt = [r*np.cos(theta), r*np.sin(theta)]
        theta += d
        pts.append(pt)
    return pts

def mk_obstacle(vertices):
    return vertices[::-1]

def mk_bounding_box(poly):
    vs = poly.vertex_list_per_poly[0]
    xs = np.sort([v[0] for i, v in vs])
    ys = np.sort([v[1] for i, v in vs])
    min_x = xs[0]
    max_x = xs[-1]
    min_y = ys[0]
    max_y = ys[-1]
    bb_verts = np.array([(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y)])
    bb = Simple_Polygon("bb"+poly.name, bb_verts)
    return min_x, max_x, min_y, max_y, bb

def uniform_sample_from_poly(poly, n):
    min_x, max_x, min_y, max_y, bb = mk_bounding_box(poly)
    samples = [[0,0]]*n
    for i in range(n):
        sample = [uniform(min_x, max_x), uniform(min_y, max_y)]
        while not IsInPoly(sample, poly):
            sample = uniform(min_x, max_x), uniform(min_y, max_y)
        samples[i] = sample
    return samples

def uniform_sample_along_circle(poly, n, r):
    samples = [[0,0]]*n
    for i in range(n):
        sample = r*normalize([uniform(-1., 1.), uniform(-1.,1.)])
        while not IsInPoly(sample, poly):
            sample = r*normalize([uniform(-1., 1.), uniform(-1.,1.)])
        samples[i] = sample
    return samples

# n: number of sides of regular polygon
# l: side length
# m: bounce trajectory hits every mth edge
# TODO: add feasibility check for theta and clockwise functionality
def regpoly_fp(n, r, m, th):
    theta = np.pi/2 - th
    c = np.cos(theta)/np.cos(theta - np.pi*(n-2*m)/n)
    l = 2*r*np.sin(np.pi/n)
    Anum = l*np.sin(np.pi*(m+1)/n)*np.sin(m*np.pi/n)
    Aden = np.sin(np.pi/n)*np.sin(np.pi*(n-2*m)/n)
    A = Anum/Aden
    xfp = (l - A*(1-c))/(1+c)
    return xfp

# returns polygon formed by collision points of limit cycle
# in regular polygon with n edges
# edge length l
# trajectory hits every mth edge of polygon
def get_limit_cycle_poly(n, r, m, th):
    poly_verts = mk_regpoly(n, r)
    edges = zip(poly_verts, (poly_verts[1:]+poly_verts[0]))
    xfp = regpoly_fp(n,r,m,th)
    cycle_points = []
    for p1, p2 in edges:
        v = normalize(np.array(p2) - np.array(p1))
        pt = np.array(p1) + v*xfp
        cycle_points.append(pt)
    return cycle_points



# Magnetic flow field generation
# all the ugly packing/unpacking is to make matplotlib stuff work
# ------------------------------

class Wire():
    def __init__(self, xy, dir):
        self.xy = xy
        self.dir = dir

    def force_at(self, x, y):
        xself, yself = self.xy
        normal = np.array([x-xself, y-yself])
        # current flowing through wire creates mag field that drops off as 1/r
        field_strength = np.true_divide(1.0, np.linalg.norm(normal, axis=0))
        field = []
        if self.dir == "CW":
            field = field_strength*normalize(rotate_vector(normal, 3.*np.pi/2.))
        elif self.dir == "CCW":
            field = field_strength*normalize(rotate_vector(normal, np.pi/2.))
        else:
            field = np.array([0.*x, 0.*y])
        return field[0], field[1]

# magnetic fields follow superposition principle
def force_from_wires(wires, xy):
    x, y = xy
    force = np.array([0.0, 0.0])
    for w in wires:
        fx, fy = w.force_at(x, y)
        force += np.array([fx, fy])
    return force

# running policies
# ----------------

# hardcoded: four wires, each wire is either CW, CCW, or X
# let wires be arranged as:
    # 0  1
    # 2  3
# let CW := 0
#     CCW := 1
#     X := 2
#
# encode with base-3 number system

wire_to_state = {"CW":0, "CCW":1, "X":2}
state_to_wire = {0:"CW", 1:"CCW", 2:"X"}

def encode_policy(wires):
    [w0, w1, w2, w3] = wires # TODO use pycontracts
    policy = 0
    for i,w in enumerate(wires):
        policy += wire_to_state[w]*(3**i)
    return policy

def decode_policy(policy):
    states = ["", "", "", ""]
    for i in range(4):
        mod_policy = policy % 3
        states[i] = mod_policy
        policy = policy // 3
    str_states = [state_to_wire[s] for s in states]
    return str_states

def encodeJointState(states):
    X = 5 # five options for state
    joint_state = 0
    N = len(states)-1
    for s in states:
        joint_state += s*(X**N)
        N = N - 1
    return joint_state
