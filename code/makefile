#CCC = icc
CCC = g++

all: sveta svtreeml svtreemu

sveta: sveta.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	#  Add -lintlc when using intel complier
	$(CCC) -std=gnu++11 sveta.cpp -o sveta -L/usr/local/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I/usr/local/include -I./gzstream  -I./lbfgsb

svtreeml: svtreeml.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) -std=gnu++11 svtreeml.cpp -o svtreeml -L/usr/local/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I. -I/usr/local/include -I./gzstream -I./lbfgsb

svtreemu: svtreemu.cpp
	cd gzstream/ && make
	cd lbfgsb/ && cmake ./ && make
	$(CCC) -std=gnu++11 svtreemu.cpp -o svtreemu -L/usr/local/lib/ -lboost_program_options -lgsl -lgslcblas -L./lbfgsb -llbfgsb -L./gzstream -lgzstream -lz -I. -I/usr/local/include -I./gzstream -I./lbfgsb


#lib:
#	$(CCC) -shared -fPIC sveta.cpp -o libsveta.so -L/usr/local/lib/ -lgsl -L./gzstream -lgzstream -I/usr/local/include

clean:
	rm -f sveta
	rm -f svtreeml
	rm -f svtreemu
	rm -f libsveta.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o
	rm -f lbfgsb/liblbfgsb.a
	rm -rf lbfgsb/CMakeFiles
	rm -f lbfgsb/Makefile
	rm -f lbfgsb/cmake_install.cmake
	rm -f lbfgsb/CMakeCache.txt
