#-------------------------------------------------------------------------------
# executable name
TARGET := libpspm

# files
SRCFILES  :=  $(wildcard src/*.cpp) 
TENSORFILES := $(wildcard tensorlib/include/*.h)
HEADERS := $(wildcard tensorlib/include/*.h) $(wildcard src/*.tpp) $(wildcard include/*.h) $(wildcard tests/*.h) 
# ------------------------------------------------------------------------------

CXX = g++
AR = ar

# paths
#CUDA_INSTALL_PATH ?= /usr/local/cuda#-5.0

# include and lib dirs (esp for cuda)
INC_PATH := -I./include -I./tensorlib/include -I/usr/include/eigen3
LIB_PATH := -L./lib 

PROFILING_FLAGS = -g -pg

# flags
CPPFLAGS = -O2 -std=c++17 -fPIC -Wall -Wextra -pedantic
LDFLAGS =  

CPPFLAGS += $(PROFILING_FLAGS)
LDFLAGS += $(PROFILING_FLAGS)

#CPPFLAGS +=   \
#-pedantic-errors -Wcast-align \
#-Wcast-qual -Wconversion \
#-Wdisabled-optimization \
#-Wfloat-equal -Wformat=2 \
#-Wformat-nonliteral -Wformat-security  \
#-Wformat-y2k \
#-Wimplicit  -Wimport  -Winit-self  -Winline \
#-Winvalid-pch   \
#-Wlong-long \
#-Wmissing-field-initializers -Wmissing-format-attribute   \
#-Wmissing-include-dirs -Wmissing-noreturn \
#-Wpacked  -Wpadded -Wpointer-arith \
#-Wredundant-decls \
#-Wshadow -Wstack-protector \
#-Wstrict-aliasing=2 -Wswitch-default \
#-Wswitch-enum \
#-Wunreachable-code -Wunused \
#-Wunused-parameter \
#-Wvariadic-macros \
#-Wwrite-strings 

CPPFLAGS += -Wno-unused-variable -Wno-sign-compare \
-Wno-unused-parameter
# -Wno-unused-but-set-variable

# libs
LIBS = 	# -ltbb #-lgsl -lgslcblas 	# additional libs
#LIBS = -lcudart 			# cuda libs 		

# files
TENSOROBJECTS := $(patsubst tensorlib/include/%.h, tensorlib/build/%.o, $(TENSORFILES))
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(SRCFILES))


all: dir $(TARGET)	

dir: 
	mkdir -p lib build tests/build tensorlib.build

$(TARGET):  $(OBJECTS) 
	$(AR) rcs lib/$(TARGET).a $(OBJECTS) $(LIBS) 
#	g++ $(LDFLAGS) -o $(TARGET) $(LIB_PATH) $(OBJECTS) $(LIBS) 

$(OBJECTS): build/%.o : src/%.cpp $(HEADERS)
	$(CXX) -c $(CPPFLAGS) $(INC_PATH) $< -o $@ 

# $(TENSOROBJECTS): tensorlib/build/%.o : tensorlib/include/%.h $(HEADERS)
# 	$(CXX) -c $(CPPFLAGS) $(INC_PATH) $< -o $@ 

clean:
	rm -f $(TARGET) build/*.o lib/*.a src/*.o
	
re: clean all

superclean: clean testclean democlean


## TESTING SUITE ##

TEST_FILES = $(wildcard tests/*.cpp) 
TEST_OBJECTS = $(patsubst tests/%.cpp, tests/%.o, $(TEST_FILES))
TEST_TARGETS = $(patsubst tests/%.cpp, tests/%.test, $(TEST_FILES))
TEST_RUNS = $(patsubst tests/%.cpp, tests/%.run, $(TEST_FILES))
ADD_OBJECTS = 

check: dir $(TARGET) compile_tests clean_log run_tests # plant_demo_test

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
	rm -f tests/*.o tests/*.test lib/*.a lib/*.so tests/1
	
dirclean:
	rm -f *.txt *.out tests/*.txt tests/*.out

recheck: testclean check

.PHONY: $(TEST_RUNS) run_tests clean testclean
# ------------------------------------------------------------------------------

demos: $(TARGET)
	cd demo/Daphnia_model && $(MAKE) FILE=fmu_equil.cpp && ./fmu_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=ifmu_equil.cpp && ./ifmu_equil.exec
#	cd demo/Daphnia_model && $(MAKE) FILE=ifmu2_equil.cpp && ./ifmu2_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=ebt_equil.cpp && ./ebt_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=iebt_equil.cpp && ./iebt_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=cm_equil.cpp && ./cm_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=icm_equil.cpp && ./icm_equil.exec
	cd demo/Daphnia_model && $(MAKE) FILE=abm_equil.cpp && ./abm_equil.exec
	cd demo/RED_model && $(MAKE) FILE=fmu_equil.cpp && ./fmu_equil.exec
	cd demo/RED_model && $(MAKE) FILE=ifmu_equil.cpp && ./ifmu_equil.exec
#	cd demo/RED_model && $(MAKE) FILE=ifmu2_equil.cpp && ./ifmu2_equil.exec
	cd demo/RED_model && $(MAKE) FILE=ebt_equil.cpp && ./ebt_equil.exec
	cd demo/RED_model && $(MAKE) FILE=iebt_equil.cpp && ./iebt_equil.exec
	cd demo/RED_model && $(MAKE) FILE=cm_equil.cpp && ./cm_equil.exec
	cd demo/RED_model && $(MAKE) FILE=icm_equil.cpp && ./icm_equil.exec
	cd demo/RED_model && $(MAKE) FILE=abm_equil.cpp && ./abm_equil.exec

democlean:
	rm -f demo/Daphnia_model/*.o demo/Daphnia_model/*.exec
	rm -f demo/RED_model/*.o     demo/RED_model/*.exec
	rm -f demo/Plant_model/*.o   demo/Plant_model/*.exec
	rm -f demo/Plant_model/src/*.o	

plant_demo_noFeedback: $(TARGET)
	cd demo/Plant_model && $(MAKE) FILE=plant_fmu_1spp.cpp && ./plant_fmu_1spp.exec && mkdir -p outputs/fmu && mv *.txt outputs/fmu
	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu_1spp.cpp && ./plant_ifmu_1spp.exec  && mkdir -p outputs/ifmu && mv *.txt outputs/ifmu
#	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu2_1spp.cpp && ./plant_ifmu2_1spp.exec  && mkdir -p outputs/ifmu2 && mv *.txt outputs/ifmu2
	cd demo/Plant_model && $(MAKE) FILE=plant_ebt_1spp.cpp && ./plant_ebt_1spp.exec  && mkdir -p outputs/ebt && mv *.txt outputs/ebt
	cd demo/Plant_model && $(MAKE) FILE=plant_iebt_1spp.cpp && ./plant_iebt_1spp.exec  && mkdir -p outputs/iebt && mv *.txt outputs/iebt
	cd demo/Plant_model && $(MAKE) FILE=plant_cm_1spp.cpp && ./plant_cm_1spp.exec  && mkdir -p outputs/cm && mv *.txt outputs/cm
	cd demo/Plant_model && $(MAKE) FILE=plant_icm_1spp.cpp && ./plant_icm_1spp.exec  && mkdir -p outputs/icm && mv *.txt outputs/icm
	cd demo/Plant_model && $(MAKE) FILE=plant_abm_1spp.cpp && ./plant_abm_1spp.exec  && mkdir -p outputs/abm && mv *.txt outputs/abm

plant_demo_test: $(TARGET)
	cd demo/Plant_model && \
	$(MAKE) FILE=plant_cm_3spp_fixedmode.cpp && \
	./plant_cm_3spp_fixedmode.exec && \
	printf "%b" "\033[0;32m[PASS]\033[m" ": plant_cm_3spp_fixedmode.cpp \n"  || \
	printf "%b" "\033[1;31m[FAIL]\033[m" ": plant_cm_3spp_fixedmode.cpp \n" && \
	Rscript plant_cm_3spp_analysis.R && \
	rm *.txt


plant_demo_withFeedback: $(TARGET)
	cd demo/Plant_model && $(MAKE) FILE=plant_fmu_1spp.cpp && ./plant_fmu_1spp.exec -1 405.32 && mkdir -p outputs/fmu_f3 && mv *.txt outputs/fmu_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu_1spp.cpp && ./plant_ifmu_1spp.exec -1 405.32 && mkdir -p outputs/ifmu_f3 && mv *.txt outputs/ifmu_f3
#	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu2_1spp.cpp && ./plant_ifmu2_1spp.exec -1 305.32 && mkdir -p outputs/ifmu2_f3 && mv *.txt outputs/ifmu2_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_ebt_1spp.cpp && ./plant_ebt_1spp.exec -1 405.32 1.0 && mkdir -p outputs/ebt_f3 && mv *.txt outputs/ebt_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_iebt_1spp.cpp && ./plant_iebt_1spp.exec -1 405.32 0.5 && mkdir -p outputs/iebt_f3 && mv *.txt outputs/iebt_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_cm_1spp.cpp && ./plant_cm_1spp.exec -1 405.32 2.0 && mkdir -p outputs/cm_f3 && mv *.txt outputs/cm_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_icm_1spp.cpp && ./plant_icm_1spp.exec -1 405.32 2.0 && mkdir -p outputs/icm_f3 && mv *.txt outputs/icm_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_abm_1spp.cpp && ./plant_abm_1spp.exec -1 105.32 50 && mkdir -p outputs/abm_f3 && mv *.txt outputs/abm_f3

plant_demo_withFeedback_IC:
#	cd demo/Plant_model && $(MAKE) FILE=plant_fmu_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_fmu_1spp.exec -1 205.32 && mkdir -p outputs/fmu_f3 && mv *.txt outputs/fmu_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_ifmu_1spp.exec -1 205.32 && mkdir -p outputs/ifmu_f3_ic && mv *.txt outputs/ifmu_f3_ic
#	cd demo/Plant_model && $(MAKE) FILE=plant_ifmu2_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_ifmu2_1spp.exec -1 205.32 && mkdir -p outputs/ifmu2_f3 && mv *.txt outputs/ifmu2_f3
#	cd demo/Plant_model && $(MAKE) FILE=plant_ebt_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_ebt_1spp.exec -1 205.32 1.0 && mkdir -p outputs/ebt_f3 && mv *.txt outputs/ebt_f3
	cd demo/Plant_model && $(MAKE) FILE=plant_iebt_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_iebt_1spp.exec -1 205.32 0.5 && mkdir -p outputs/iebt_f3_ic && mv *.txt outputs/iebt_f3_ic
#	cd demo/Plant_model && $(MAKE) FILE=plant_cm_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_cm_1spp.exec -1 205.32 2.0 && mkdir -p outputs/cm_f3 && mv *.txt outputs/cm_f3
#	cd demo/Plant_model && $(MAKE) FILE=plant_icm_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_icm_1spp.exec -1 205.32 2.0 && mkdir -p outputs/icm_f3_ic && mv *.txt outputs/icm_f3_ic
#	cd demo/Plant_model && $(MAKE) FILE=plant_abm_1spp.cpp INIT_DIST=-DUSE_INIT_DIST && ./plant_abm_1spp.exec -1 205.32 50 && mkdir -p outputs/abm_f3 && mv *.txt outputs/abm_f3


website:
	R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); pkgdown::clean_site(); pkgdown::init_site(); pkgdown::build_home(); pkgdown::build_articles(); pkgdown::build_tutorials(); pkgdown::build_news()"

api:
	doxygen doxygen/Doxyfile




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
	
