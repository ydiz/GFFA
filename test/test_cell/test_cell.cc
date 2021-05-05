#include "../GFFA.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // Coordinate latt_size   = GridDefaultLatt();
  Coordinate latt_size   = Coordinate({6,6,6,6});
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();

  Coordinate cell_size = Coordinate({4,4,4,4});

  GridCartesian             Grid(latt_size, simd_layout, mpi_layout);
  GridCartesian             cell_Grid(cell_size, simd_layout, mpi_layout);

  LatticeGaugeField P(&Grid); P = 1.0;


  // Test 1. applying mask to P
  LatticeLorentzScalar mask = get_mask(P.Grid(), cell_size);
  // std::cout << mask << std::endl;
  P = P * mask;
  // std::cout << P << std::endl;
  
  // Test 2. inserting a smaller cell into a larger lattice
  LatticeGaugeField cell_P(&cell_Grid); cell_P = 2.0;
  localCopyRegion(cell_P, P, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

  std::cout << P << std::endl;



  Grid_finalize();
}


