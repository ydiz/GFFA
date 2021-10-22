#include "../GFFA.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // Coordinate latt_size   = GridDefaultLatt();
  Coordinate latt_size   = Coordinate({16,16,16,16});
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  Coordinate mpi_size  = GridDefaultMpi();
  GridCartesian             grid(latt_size, simd_layout, mpi_size);

  Coordinate cell_size   = Coordinate({6,6,6,6});
  Coordinate cell_grid_size(4);
  for(int i=0; i<4; ++i) cell_grid_size[i] = mpi_size[i] * cell_size[i];
  GridCartesian cell_grid(cell_grid_size, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_size);

  
  LatticeGaugeField U(&grid);
  FieldMetaData header;
  NerscIO::readConfiguration(U, header, "ckpoint_lat.500");

  LatticeGaugeField U_cell(&cell_grid); U_cell = Zero();
  localCopyRegion(U, U_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

  std::cout << "print U" << std::endl;
  print_grid_field_site(U, {0,0,0,0});
  std::cout << "print U_cell" << std::endl;
  print_grid_field_site(U_cell, {0,0,0,0});


  Grid_finalize();
}


