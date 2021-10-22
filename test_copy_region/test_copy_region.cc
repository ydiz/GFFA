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

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > coor(&grid);
  for(int mu=0;mu<4;mu++){
    LatticeComplex tmp(&grid);
    LatticeCoordinate(tmp, mu);
    pokeLorentz(coor, tmp, mu);
  }

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > cell_coor(&cell_grid);
  localCopyRegion(coor, cell_coor, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

  // std::cout << cell_coor << std::endl;

  GridRedBlackCartesian rbGrid(cell_coor.Grid());
  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > cell_coor_oe(&rbGrid);
  pickCheckerboard(Even, cell_coor_oe, cell_coor);
  print_half_field(cell_coor_oe);

  Grid_finalize();
}


