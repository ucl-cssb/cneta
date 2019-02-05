# CCC = clang++ -O0
CCC = g++

# Change the value of system and directory accordingly
SYSTEM     = x86-64_osx
LIBFORMAT  = static_pic
CPLEXDIR      = /Applications/CPLEX_Studio128/cplex
CONCERTDIR    = /Applications/CPLEX_Studio128/concert

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -stdlib=libc++

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
CCLNDIRS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)

CCLNFLAGS = -lconcert -lilocplex -lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit -stdlib=libc++

GSLFLAGS = -lgsl -lgslcblas

mcomp_ilp: mcomp_ilp.o
	$(CCC) $(GSLFLAGS) $(CCFLAGS) $(CCLNDIRS) -o mcomp_ilp mcomp_ilp.o $(CCLNFLAGS)
mcomp_ilp.o: mcomp_ilp.cpp
	$(CCC) -c $(CCFLAGS) mcomp_ilp.cpp -o mcomp_ilp.o

clean:
	rm  mcomp_ilp mcomp_ilp.o 
