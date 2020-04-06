import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set(rc={'figure.figsize':(11,4)})
from copy import copy
from time import sleep

from utilities import *
from configuration import *

def calc_centroid(points, k):
    centroid = np.sum(points)/k
    return centroid

def find_densities(fname, run):
    densities = []
    with open(fname,'r') as f:
        t = 0
        for line in f:
            it = iter([float(x) for x in line.split()])
            coords = np.array(list(zip(it, it))) # grab consecutive pairs. python black magic
            IN, OUT = 0, 0
            old_centroid = np.array([0., 0.])
            k = coords.shape[0]
            centroid = calc_centroid(coords, k)

            ds = np.empty(k)
            for i,pt in enumerate(coords):
                d = np.linalg.norm(np.array(pt) - centroid)
                ds[i] = d
            densities.append([run, t,np.average(ds)])
            t+=1

    return pd.DataFrame(densities, columns=["Run","T","D"])


baselines = ["data/oct_baseline.xyz",
             "data/oct_baseline_2.xyz",
             "data/oct_baseline_3.xyz",
             "data/oct_baseline_4.xyz",
             "data/oct_baseline_5.xyz"]
tens =  ["data/oct_10percent.xyz",
         "data/oct_10percent_2.xyz",
         "data/oct_10percent_3.xyz",
         "data/oct_10percent_4.xyz",
         "data/oct_10percent_5.xyz"]
twenties =  ["data/oct_20percent.xyz",
         "data/oct_20percent_2.xyz",
         "data/oct_20percent_3.xyz",
         "data/oct_20percent_4.xyz",
         "data/oct_20percent_5.xyz"]
thirties =  ["data/oct_30percent.xyz",
         "data/oct_30percent_2.xyz",
         "data/oct_30percent_3.xyz",
         "data/oct_30percent_4.xyz",
         "data/oct_30percent_5.xyz"]


if __name__ == '__main__':

    print("Analyzing baseline data...")
    baseline_densities = [find_densities(f,i) for i,f in enumerate(baselines)]
    bds = pd.concat(baseline_densities)
    print(bds)
    #bds_lumped = bds.agg(['mean','std'], axis="columns")
    #print(bds_lumped)
    sns.lineplot(x=bds["T"], y=bds["D"],hue=bds["Run"], estimator=None)
    plt.show()
    #baseline_densities = [find_densities(f) for f in baselines]
    #print("Analyzing 10% data...")
    #ten_percent_densities = [find_densities(f) for f in tens]
    #print("Analyzing 20% data...")
    #twenty_percent_densities = [find_densities(f) for f in twenties]
    #print("Analyzing 30% data...")
    #thirty_percent_densities = [find_densities(f) for f in thirties]
    #print("Rendering Plot...")
    #t = range(0, len(baseline_densities))
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plt.plot(t, baseline_densities, '-k', linewidth=2, label='0%')
    #plt.plot(t, ten_percent_densities, '-r', linewidth=2, label='10%')
    #plt.plot(t, twenty_percent_densities, '-g', linewidth=2, label='20%')
    #plt.plot(t, thirty_percent_densities, '-b', linewidth=2, label='30%')
    #ax.legend()
    #plt.xlabel("Timestep")
    #plt.ylabel("Average Distance to Centroid")
    #plt.savefig("density.pdf", bbox='tight')

## computes theta resulting in counterclockwise trajectory in middle of
## admissable range
## TODO: add clockwise functionality
#def theta_stable(n, l, m):
#    phi = (n-2)*np.pi/n
#    phi_m = np.pi*(n-2*m)/n
#    phi_mless1 = np.pi*(n-2*(m-1))/n
#    theta = (phi_m + phi_mless1)/4.
#    return theta
#
#
#
## area of regular polygon
## limit cycle of regular polygon will be regular polyon
## TODO: general area calculator (from triangulation)
#def poly_area(n, r):
#    area = (r**2)*n*np.sin(2*np.pi/n)/2
#    return area
