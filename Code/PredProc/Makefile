CXX = g++
OPTFLAGS = -O3 -fopenmp
CXXFLAGS = -I/opt/OpenBLAS/include 
LIBS = -I/opt/OpenBLAS/include -I/usr/local/lib/R/include -L /opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran -L/usr/local/lib -lRmath -L/usr/local/lib/R/lib -lR

gpp: gpp.o kvpar.o
	$(CXX) $(OPTFLAGS) -o gpp gpp.o kvpar.o $(LIBS) 

gpp.o: gpp.cpp
	$(CXX) $(OPTFLAGS) gpp.cpp -c $(LIBS)

kvpar.o: kvpar.cpp
	$(CXX) $(OPTFLAGS) kvpar.cpp -c

clean:
	rm -rf *o gpp
