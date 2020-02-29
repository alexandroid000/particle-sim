
# vary frac in 10%, 20%, 30%

# 300 agents
AGENTS=300

# 400 time steps
T=400

for i in {1..3}
do
    ( source psim/bin/activate \
    ./run_sim.py $i & )
done

