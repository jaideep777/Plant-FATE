FILE=ebtn_equil
g++ -O3 -g -pg -I../../include -c ../../src/species.cpp ../../src/rkck45.cpp ../../src/lsoda.cpp ../../src/solver.cpp ../../src/ebt.cpp  $FILE.cpp
g++ -g -pg -o 1 species.o rkck45.o lsoda.o solver.o ebt.o $FILE.o

./1 && \
	printf "%b" "\033[0;32m[PASS]\033[m" ": $FILE \n"  || \
	printf "%b" "\033[1;31m[FAIL]\033[m" ": $FILE \n"
