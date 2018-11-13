# for hadron, add -fpermissive to CXXFLAGS

CXX = mpicxx

CXXFLAGS =  -I/home/ydzhao/cuth/install/Grid/include -I/home/ydzhao/cuth/install/openmpi-3.1.0/include -I/home/ydzhao/cuth/install/hdf5-1.8.20/include -I/home/ydzhao/cuth/install/LIME/include -I/home/ydzhao/cuth/install/fftw-3.3.8/include -fopenmp -O3 -std=c++11

LDFLAGS = -L/home/ydzhao/cuth/install/Grid/lib -L/home/ydzhao/cuth/install/openmpi-3.1.0/lib -Wl,-rpath -Wl,/home/ydzhao/cuth/install/openmpi-3.1.0/lib -L/home/ydzhao/cuth/install/hdf5-1.8.20/lib -L/home/ydzhao/cuth/install/LIME/lib -L/home/ydzhao/cuth/install/fftw-3.3.8/lib -fopenmp

LIBS = -lGrid -lmpi_cxx -lmpi -Wl,-rpath -Wl,/home/ydzhao/cuth/install/openmpi-3.1.0/lib  -lhdf5_cpp -lz -lcrypto -llime -lfftw3f -lfftw3 -lmpfr -lgmp -lstdc++ -lm -lrt -lz -lhdf5

TARGET = Test_hmc

$(TARGET): FORCE
	$(CXX) $(CXXFLAGS) -c -o $(TARGET).o $(TARGET).cc
	$(CXX) $(LDFLAGS) -o $@ $(TARGET).o $(LIBS)
    
clean:
	-rm *.o

.PHONY: FORCE
