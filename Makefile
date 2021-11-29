#-------------------------------------------------------------------------------
# executable name
TARGET := 1

# files
SRCFILES  :=  $(wildcard src/*.cpp) 
HEADERS := $(wildcard src/*.tpp) $(wildcard include/*.h) $(wildcard tests/*.h)
# ------------------------------------------------------------------------------

# paths
#CUDA_INSTALL_PATH ?= /usr/local/cuda#-5.0

# include and lib dirs (esp for cuda)
INC_PATH :=  -I./include #-I./CppNumericalSolvers-1.0.0
INC_PATH += -I/home/jaideep/codes/phydro_cpp/include -I/home/jaideep/codes/phydro_cpp/LBFGSpp/include -I/usr/include/eigen3 
LIB_PATH :=  

# flags
CPPFLAGS = -O3 -g -pg -std=c++11 -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable
LDFLAGS =  -g -pg

# libs
LIBS = 	 -lgsl -lgslcblas 	# additional libs
#LIBS = -lcudart 			# cuda libs 		

# files
OBJECTS = $(patsubst src/%.cpp, build/%.o, $(SRCFILES))


all: dir $(TARGET)	

dir: 
	mkdir -p lib build tests/build

$(TARGET): $(OBJECTS) 
	g++ $(LDFLAGS) -o $(TARGET) $(LIB_PATH) $(OBJECTS) $(LIBS) 

$(OBJECTS): build/%.o : src/%.cpp $(HEADERS)
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@ 

clean:
	rm -f $(TARGET) build/*.o log.txt gmon.out 
	
re: clean all

superclean: clean testclean


## TESTING SUITE ##

TEST_FILES = $(wildcard tests/*.cpp) 
TEST_OBJECTS = $(patsubst tests/%.cpp, tests/%.o, $(TEST_FILES))
TEST_TARGETS = $(patsubst tests/%.cpp, tests/%.test, $(TEST_FILES))
TEST_RUNS = $(patsubst tests/%.cpp, tests/%.run, $(TEST_FILES))
ADD_OBJECTS = 

check: $(OBJECTS) compile_tests clean_log run_tests

compile_tests: $(TEST_TARGETS)
	
clean_log:
	@rm -f log.txt

run_tests: $(TEST_RUNS) 
	
$(TEST_RUNS): tests/%.run : tests/%.test	
	@echo "~~~~~~~~~~~~~~~ $< ~~~~~~~~~~~~~~~~" >> log.txt
	@./$< >> log.txt && \
		printf "%b" "\033[0;32m[PASS]\033[m" ": $* \n"  || \
		printf "%b" "\033[1;31m[FAIL]\033[m" ": $* \n"

$(TEST_OBJECTS): tests/%.o : tests/%.cpp $(HEADERS) 
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@

$(TEST_TARGETS): tests/%.test : tests/%.o
	g++ $(LDFLAGS) -o $@ $(LIB_PATH) $(OBJECTS) $(ADD_OBJECTS) $< $(LIBS) 

testclean: 
	rm -f tests/*.o tests/*.test

recheck: testclean check

.PHONY: $(TEST_RUNS) run_tests clean testclean
# ------------------------------------------------------------------------------





#-gencode=arch=compute_10,code=\"sm_10,compute_10\"  -gencode=arch=compute_20,code=\"sm_20,compute_20\"  -gencode=arch=compute_30,code=\"sm_30,compute_30\" 

#-W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function 
#-m64 -fno-strict-aliasing 
#-I. -I/usr/local/cuda/include -I../../common/inc -I../../../shared//inc 
#-DUNIX -O2


#g++ -fPIC   -m64 -o ../../bin/linux/release/swarming_chasing_predator obj/x86_64/release/genmtrand.cpp.o  obj/x86_64/release/simpleGL.cu.o  -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart   
#-lGL -lGLU -lX11 -lXi -lXmu -lGLEW_x86_64 -L/usr/X11R6/lib64 -lGLEW_x86_64 -L/usr/X11R6/lib64 -lglut 
#-L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart 
#-L/usr/lib -lgsl -lgslcblas 
#-lcutil_x86_64  -lshrutil_x86_64 




#CXXWARN_FLAGS := \
#	-W -Wall \
#	-Wimplicit \
#	-Wswitch \
#	-Wformat \
#	-Wchar-subscripts \
#	-Wparentheses \
#	-Wmultichar \
#	-Wtrigraphs \
#	-Wpointer-arith \
#	-Wcast-align \
#	-Wreturn-type \
#	-Wno-unused-function \
#	$(SPACE)

#CWARN_FLAGS := $(CXXWARN_FLAGS) \
#	-Wstrict-prototypes \
#	-Wmissing-prototypes \
#	-Wmissing-declarations \
#	-Wnested-externs \
#	-Wmain \
#	
	
#HEADERS  := $(wildcard *.h)
	
