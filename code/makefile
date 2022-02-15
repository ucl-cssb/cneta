# CCC = icc		# Using intel complier
CCC = g++
# CCC = clang++   # when compiling on Mac
BOOST = /usr/local
# BOOST = /share/apps/boost-1.71.0
# GSL = /share/apps/gsl-2.5/
GSL = /usr/local
omp = -fopenmp
# omp =
FLAG = -O3 -std=gnu++11
# FLAG =
# FLAG = -O3 -lintlc 	#  Add -lintlc when using intel complier

all: cnets cnetml cnetmcmc 

cnets: cnets.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) cnets.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp evo_tree.cpp tree_op.cpp genome.cpp model.cpp parse_cn.cpp stats.cpp -o cnets -L$(BOOST)/lib/ -lboost_program_options -L$(GSL)/lib/ -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I$(BOOST)/include -I$(GSL)/include -I./gzstream  -I./lbfgsb

cnetml: cnetml.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) $(omp) cnetml.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp stats.cpp evo_tree.cpp tree_op.cpp model.cpp likelihood.cpp nni.cpp optimization.cpp parse_cn.cpp state.cpp -o cnetml -L$(BOOST)/lib/ -lboost_filesystem -lboost_system -lboost_program_options -L$(GSL)/lib/ -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I./ -I$(BOOST)/include -I$(GSL)/include -I./gzstream -I./lbfgsb

cnetmcmc: cnetmcmc.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) cnetmcmc.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp stats.cpp evo_tree.cpp tree_op.cpp model.cpp likelihood.cpp nni.cpp optimization.cpp parse_cn.cpp -o cnetmcmc -L$(BOOST)/lib/ -lboost_program_options -L$(GSL)/lib/ -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz  -I./ -I$(BOOST)/include -I$(GSL)/include -I./gzstream -I./lbfgsb

# test: test.cpp
#	cd gzstream/ && make
#	cd lbfgsb/ && cmake ./ && make
#	$(CCC) $(FLAG) $(omp) test.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp stats.cpp evo_tree.cpp tree_op.cpp model.cpp likelihood.cpp nni.cpp optimization.cpp parse_cn.cpp state.cpp -o test -L$(BOOST)/lib/ -lboost_filesystem -lboost_system -lboost_program_options -L$(GSL)/lib/ -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I./ -I$(BOOST)/include -I$(GSL)/include -I./gzstream -I./lbfgsb

#lib:
#	$(CCC) -shared -fPIC cnets.cpp -o libcnets.so -L$(BOOST)/lib/ -lgsl -L./gzstream -lgzstream -I$(BOOST)/include

clean:
	rm -f cnets
	rm -f cnetml
	rm -f cnetmcmc
	rm -f libcnets.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o
	rm -f lbfgsb/liblbfgsb.a
	rm -rf lbfgsb/CMakeFiles
	rm -f lbfgsb/Makefile
	rm -f lbfgsb/cmake_install.cmake
	rm -f lbfgsb/CMakeCache.txt
