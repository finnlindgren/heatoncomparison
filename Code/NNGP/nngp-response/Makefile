CXX = g++
OPTFLAGS = -O3 -fopenmp
CXXFLAGS = -I/opt/OpenBLAS/include 
LIBS =  -I/home/finleya2/local/include/ -L/home/finleya2/local/lib/ -lopenblas -lpthread -lgfortran -lRmath -L/home/finleya2/local/lib/R/lib -lR -I/home/finleya2/local/lib/R/include/

nngp: nngp.o kvpar.o
	$(CXX) $(OPTFLAGS) -o nngp nngp.o kvpar.o $(LIBS) 

nngp.o: nngp.cpp
	$(CXX) $(OPTFLAGS) nngp.cpp -c $(LIBS)

kvpar.o: libs/kvpar.cpp
	$(CXX) $(OPTFLAGS) libs/kvpar.cpp -c

clean:
	rm -rf *o nngp
