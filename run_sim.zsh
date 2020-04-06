
# vary frac in 10%, 20%, 30%

# 300 agents
AGENTS=300

# 400 time steps
T=400

for i in {0..3}
do
    echo "running" $i
    source psim/bin/activate
    ./run_sim.py $i &
done

wait

mv octagon_N300_T400_F0.0_typeA.xyz data/oct_baseline_5.xyz
mv octagon_N300_T400_F0.1_typeA.xyz data/oct_10percent_5.xyz
mv octagon_N300_T400_F0.2_typeA.xyz data/oct_20percent_5.xyz
mv octagon_N300_T400_F0.3_typeA.xyz data/oct_30percent_5.xyz
