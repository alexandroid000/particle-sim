import numpy as np
import matplotlib.pyplot as plt
from utilities import *
from configuration import *
from copy import copy
from time import sleep
from simple_polygon import Simple_Polygon
from helper.shoot_ray_helper import IsInPoly
from utilities import mk_regpoly
from environments import octagon

# computes theta resulting in counterclockwise trajectory in middle of
# admissable range
# TODO: add clockwise functionality
def theta_stable(n, l, m):
    phi = (n-2)*np.pi/n
    phi_m = np.pi*(n-2*m)/n
    phi_mless1 = np.pi*(n-2*(m-1))/n
    theta = (phi_m + phi_mless1)/4.
    return theta

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

# area of regular polygon
# limit cycle of regular polygon will be regular polyon
# TODO: general area calculator (from triangulation)
def poly_area(n, r):
    area = (r**2)*n*np.sin(2*np.pi/n)/2
    return area

def find_densities(fname, frac, c_poly, area_c, area_outside):
    densities = []
    with open(fname,'r') as f:
        for line in f:
            it = iter([float(x) for x in line.split()])
            coords = zip(it, it) # grab consecutive pairs. python black magic
            IN, OUT = 0, 0
            for pt in coords:
                if IsInPoly(pt, c_poly):
                    IN += 1
                else: # assuming no escape polygon...
                    OUT += 1
            in_density = IN/(300*frac)
            densities.append(in_density)
    return densities


if __name__ == '__main__':

    theta = 0.2
    cycle = get_limit_cycle_poly(8, 1., 1, theta)
    c_poly = Simple_Polygon("cycle", cycle)
    l_cycle = np.linalg.norm(np.array(cycle[0])-np.array(cycle[1]))
    r = 1./(2.*np.sin(np.pi/8.))
    env = octagon(r)

    area_e = poly_area(8., 1.)
    area_c = poly_area(8., l_cycle)
    area_outside = area_e - area_c

    print("Analyzing baseline data...")
    baseline_densities = find_densities("data/octagon_N300_T300_F0.0_typeA.xyz", 1., c_poly, area_c, area_outside)
    print("Analyzing 10% data...")
    ten_percent_densities = find_densities("data/octagon_N300_T300_F0.1_typeA.xyz", 0.9, c_poly, area_c, area_outside)
    print("Analyzing 20% data...")
    twenty_percent_densities = find_densities("data/octagon_N300_T300_F0.2_typeA.xyz", 0.8, c_poly, area_c, area_outside)
    print("Rendering Plot...")
    t = range(0, len(baseline_densities))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(t, baseline_densities, '-g', linewidth=2)
    plt.plot(t, ten_percent_densities, '-b', linewidth=2)
    plt.plot(t, twenty_percent_densities, '-r', linewidth=2)
    plt.xlabel("Timestep")
    plt.ylabel("Fraction of Brownian Particles Inside Limit Cycle")
    plt.savefig("density.pdf", bbox='tight')


#nx, ny = 10, 10 
#y = np.linspace(-L, L, ny)
#X, Y = np.meshgrid(x, y)
#
#
#
#
#
#
#
#
#
#R = 0.85*L
#
#ax.set_xlim(-R,R)
#ax.set_ylim(-R,R)
#ax.set_aspect('equal')
#plt.savefig("policy25field.eps", bbox='tight')
#
#times = [4.27, 21.72, 119.4, 644.17, 3413.6, 18626.6, 98147.5, 522008.7]
#
#time_dat = zip(range(2,10), times)
#
#plt.clf()
#ax = fig.add_subplot(111)
