gprof tests/fmu_equil.test | ./gprof2dot.py -ws -n0.01 -e0.01 | dot -Tpng -o test_pspm_iFMU.png
