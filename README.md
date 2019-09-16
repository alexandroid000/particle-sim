#### Dependencies

- Python3 

The necessary libraries are included in requirements.txt. To automatically
install all the packages, run

```
pip install -r requirements.txt
```

### Getting Started

First, add the submodule:

```
git submodule init
git submodule update
```

Then, try running

```
python run_sim.py
```

This should write three files, a trajectory log ending in `.xyz`, a log of how
many particles were in each region ending in `.csv`, and a video of the
simulation, ending in `.mp4`.

### Contributors

