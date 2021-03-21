# CCC = icc		# Using intel complier
CCC = g++
# CCC = clang
BOOST = /usr/local
omp = -fopenmp
# omp =
FLAG = -O3
# FLAG =

all: sveta svtreeml svtreemu svtreemcmc

sveta: sveta.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	#  Add -lintlc when using intel complier
	$(CCC) $(FLAG) -std=gnu++11 sveta.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp -o sveta -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I$(BOOST)/include -I./gzstream  -I./lbfgsb

svtreeml: svtreeml.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) -std=gnu++11 $(omp) svtreeml.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp -o svtreeml -L$(BOOST)/lib/ -lboost_filesystem -lboost_system -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I./ -I$(BOOST)/include -I./gzstream -I./lbfgsb

svtreemu: svtreemu.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) -std=gnu++11 svtreemu.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp -o svtreemu -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I./ -I$(BOOST)/include -I./gzstream -I./lbfgsb

svtreemcmc: svtreemcmc.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) $(FLAG) -std=gnu++11 svtreemcmc.cpp matexp/matrix_exponential.cpp matexp/r8lib.cpp -o svtreemcmc -L$(BOOST)/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz  -I./ -I$(BOOST)/include -I./gzstream -I./lbfgsb

#lib:
#	$(CCC) -shared -fPIC sveta.cpp -o libsveta.so -L$(BOOST)/lib/ -lgsl -L./gzstream -lgzstream -I$(BOOST)/include

clean:
	rm -f sveta
	rm -f svtreeml
	rm -f svtreemu
	rm -f svtreemcmc
	rm -f libsveta.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o
	rm -f lbfgsb/liblbfgsb.a
	rm -rf lbfgsb/CMakeFiles
	rm -f lbfgsb/Makefile
	rm -f lbfgsb/cmake_install.cmake
	rm -f lbfgsb/CMakeCache.txt
