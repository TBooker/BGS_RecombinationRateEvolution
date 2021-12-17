
cd simple_map_constant_s0.01
parallel "slim -d R={1} -d s=-0.01 -d REP={2} ../configs/recombinationRateEvolution.singleElement.slim" ::: 0.1 1. 10. ::: $(seq 51 75)
cd ../

python bin/parseTreesSimple.py --tree simple_map_constant_s0.01/ --output summary_simple_map_constant_s0.01.2.csv --nProc 5
exit 0

cd simple_map_constant_s0.001
parallel "slim -d R={1} -d s=-0.001 -d REP={2} ../configs/recombinationRateEvolution.singleElement.slim" ::: 0.1 1. 10. ::: $(seq 1 50)
cd ../

python bin/parseTreesSimple.py --tree simple_map_constant_s0.001/ --output summary_simple_map_constant_s0.001.csv --nProc 5


cd simple_map_constant_s0.0002
parallel "slim -d R={1} -d s=-0.0002 -d REP={2} ../configs/recombinationRateEvolution.singleElement.slim" ::: 0.1 1. 10. ::: $(seq 1 50)
cd ../

python bin/parseTreesSimple.py --tree simple_map_constant_s0.0002/ --output summary_simple_map_constant_s0.0002.csv --nProc 5
