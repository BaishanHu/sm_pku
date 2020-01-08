exe=test.exe
#CXX=icpc -I/opt/gopenmpi/include -pthread -Wl,-rpath -Wl,/opt/gopenmpi/lib -Wl,--enable-new-dtags -L/opt/gopenmpi/lib -lmpi_cxx -lmpi
#CXX=/opt/gopenmpi/bin/mpic++
#CXX=/opt/mpich-3.2/bin/mpic++
CXX=mpic++
#CXXDEF= -D USE_MPI
#IDIR=-I /opt/include/ -I/opt/boost_1_56_0/ -I/opt/gopenmpi/include/ 
IDIR=-I /opt/include/ -I/opt/boost_1_68_0/include/
LDIR=-L /opt/boost_1_68_0/lib #-L/opt/lib/
#CXXFLAGS= -O3 #-Wall
CXXFLAGS=  -O3 -fopenmp  -std=c++11 -Wno-deprecated-declarations -Wl,--rpath=/opt/boost_1_68_0/lib
#CXXFLAGS=  -O3 -qopenmp  -std=c++11 -Wl,--rpath=/home_data/hubsh/usr/lib
#CXXFLAGS=-g -std=c++11 
GSL_LIBS=-lgsl -lgslcblas 
BOOST_LIBS=-lboost_serialization -lboost_mpi
#BOOST_LIBS=-lboost_serialization -lboost_timer -lboost_system

DEPS=help.h\
base.h\
orbital.h\
twobodystate.h\
system_ho.h\
sm_system.h\
serialization.h\
lanczos_mpi.h\
sm_solver.h

OBJ= main.o help.o orbital.o twobodystate.o system_ho.o sm_system.o sm_solver.o

%.o: %.cpp $(DEPS)
	$(CXX) ${CXXDEF} -c -o $@ $< $(CXXFLAGS) $(IDIR) $(LDIR)
$(exe):$(OBJ)
	$(CXX)  -o $@ $^ $(CXXFLAGS) $(GSL_LIBS) $(BOOST_LIBS) $(IDIR) $(LDIR) 

.PHONY: clean

clean:
	rm -f *.o  $(exe)

