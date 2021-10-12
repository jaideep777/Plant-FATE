gprof tests/test_LSODA.test | ./gprof2dot.py -ws -n0.01 -e0.01 | dot -Tpng -o test_LSODA.png
