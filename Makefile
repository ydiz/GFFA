# fftw is installed in /home/ckelly/knl/install_base_mpi
CXX = mpiicpc

CXXFLAGS = -DHAVE_CONFIG_H -I. -I/home/ydzhao/knl/install/Grid/include -xmic-avx512 -I/home/ckelly/knl/install_base_mpi/include -fopenmp  -O3 -g -std=c++11

LDFLAGS = -fopenmp -O3 -g -std=c++11 -L/home/ydzhao/knl/install/Grid/lib -L/home/ckelly/knl/install_base_mpi/lib

LIBS =  -lGrid -lz -llime -lmpfr -lgmp -lmkl_rt -lstdc++ -lm

TARGET = Test_hmc

$(TARGET): FORCE
	$(CXX) $(CXXFLAGS) -c $@.cc
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

Test_SGF: FORCE
	$(CXX) $(CXXFLAGS) -c $@.cc
	$(CXX) $(LDFLAGS) -o $@ $@.o $(LIBS)

clean:
	-rm *.o

FORCE:

.PHONY: FORCE
