#include "../GFFA.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // Coordinate latt_size   = GridDefaultLatt();
  Coordinate latt_size   = Coordinate({8,8,8,8});
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();


  GridCartesian             Grid(latt_size, simd_layout, mpi_layout);

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > coor(&Grid);

  for(int mu=0;mu<4;mu++){
    LatticeComplex tmp(&Grid);
    LatticeCoordinate(tmp, mu);
    pokeLorentz(coor, tmp, mu);
  }

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > shifted(&Grid);
  int dimension = 3, shift= 1;
  Cshift_local(shifted, coor, dimension, shift);

  std::cout << shifted << std::endl;

  Grid_finalize();
}


