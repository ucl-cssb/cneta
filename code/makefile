# CCC = icc		# Using intel complier
CCC = g++
# CCC = clang++   # when compiling on Mac
BOOST = /usr/local
# omp = -fopenmp
omp =
FLAG = -O3 -std=gnu++11
# FLAG =
# FLAG = -O3 -lintlc 	#  Add -lintlc when using intel complier

all: sveta svtreeml svtreemcmc

sveta: sveta.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) sveta.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp evo_tree.cpp tree_op.cpp genome.cpp model.cpp parse_cn.cpp stats.cpp -o sveta -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I$(BOOST)/include -I./gzstream  -I./lbfgsb

svtreeml: svtreeml.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) $(omp) svtreeml.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp stats.cpp evo_tree.cpp tree_op.cpp model.cpp likelihood.cpp nni.cpp optimization.cpp parse_cn.cpp state.cpp -o svtreeml -L$(BOOST)/lib/ -lboost_filesystem -lboost_system -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I./ -I$(BOOST)/include -I./gzstream -I./lbfgsb

svtreemcmc: svtreemcmc.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) svtreemcmc.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp stats.cpp evo_tree.cpp tree_op.cpp model.cpp likelihood.cpp nni.cpp optimization.cpp parse_cn.cpp -o svtreemcmc -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz  -I./ -I$(BOOST)/include -I./gzstream -I./lbfgsb

#lib:
#	$(CCC) -shared -fPIC sveta.cpp -o libsveta.so -L$(BOOST)/lib/ -lgsl -L./gzstream -lgzstream -I$(BOOST)/include

clean:
	rm -f sveta
	rm -f svtreeml
	rm -f svtreemcmc
	rm -f libsveta.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o
	rm -f lbfgsb/liblbfgsb.a
	rm -rf lbfgsb/CMakeFiles
	rm -f lbfgsb/Makefile
	rm -f lbfgsb/cmake_install.cmake
	rm -f lbfgsb/CMakeCache.txt
