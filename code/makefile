
all: sveta svtreeml svtreemu

sveta: sveta.cpp
	cd gzstream/ && make
	g++ -std=gnu++11 sveta.cpp -o sveta -L/usr/local/lib/ -lgsl -lgslcblas -L./gzstream -lgzstream -lz -I/usr/local/include -I./gzstream

svtreeml: svtreeml.cpp
	cd gzstream/ && make
	g++ -std=gnu++11 svtreeml.cpp -o svtreeml -L/usr/local/lib/ -lboost_program_options -lgsl -lgslcblas -L./gzstream -lgzstream -lz -I. -I/usr/local/include -I./gzstream

svtreemu: svtreemu.cpp
	cd gzstream/ && make
	g++ -std=gnu++11 svtreemu.cpp -o svtreemu -L/usr/local/lib/ -lboost_program_options -lgsl -lgslcblas -L./gzstream -lgzstream -lz -I. -I/usr/local/include -I./gzstream


#./sveta ./ 5 1

#lib:
#	g++ -shared -fPIC sveta.cpp -o libsveta.so -L/usr/local/lib/ -lgsl -L./gzstream -lgzstream -I/usr/local/include


clean:
	rm -f sveta
	rm -f svtreeml
	rm -f svtreemu
	rm -f libsveta.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o
