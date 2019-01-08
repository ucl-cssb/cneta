
sveta: sveta.cpp
	cd gzstream/ && make
	g++ sveta.cpp -o sveta -L/usr/local/lib/ -lgsl -L./gzstream -lgzstream -lz -I/usr/local/include -I./gzstream

#./sveta ./ 5 1

#lib:
#	g++ -shared -fPIC sveta.cpp -o libsveta.so -L/usr/local/lib/ -lgsl -L./gzstream -lgzstream -I/usr/local/include


clean:
	rm -f sveta
	rm -f libsveta.o
	rm -f gzstream/libgzstream.a
	rm -f gzstream/gzstream.o 
