gprof tests/fmu_equil.test | ./gprof2dot.py -ws -n0.01 -e0.01 | dot -Tpng -o test_single_plant_cm.png
