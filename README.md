#### Dependencies and Installation

- Python3 
- matplotlib
- bounce-viz, another github repo for geometry calculations (included as
submodule)

The necessary libraries are included in requirements.txt. To automatically
install the dependencies, run

```
pip install -r requirements.txt
```

Then, add the submodule:

```
git submodule init
git submodule update
```

### Getting Started


To execute a simulation, run

```
python run_sim.py
```

This should write three files, a trajectory log ending in `.xyz`, a log of how
many particles were in each region ending in `.csv`, and a video of the
simulation, ending in `.mp4`.

### Changing Parameters

The changeable aspects of the simulation are stored in `configuration.py`, which
is documented in-line. This is the only file that will usually need to be
changed to run a new simulation.
