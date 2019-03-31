#! /usr/bin/sh

# export PYTHONPATH="${PYTHONPATH}":"${PWD}"

# to collect transition data, set N=100 and T=30 in the simulation config 
# (instead of T=30, can use whatever timescale you're interested in)
# start all agents in the region of interest (ie: uniformly distributed in start state)
# look in simname_regions.csv for # of agents in all states at the end of the time window

for start in {0..4}; do
    for action in {0..81}; do
    	python3 run_sim.py "$start" "$action"
    done
done
