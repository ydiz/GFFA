#include "../GFFA.h"


NAMESPACE_BEGIN(Grid);
template<class vobj> Lattice<vobj> my_Cshift_local_rb(const Lattice<vobj> &rhs, int dimension,int shift)
{
  Lattice<vobj> ret(rhs.Grid());
  ret.Checkerboard() = rhs.Grid()->CheckerBoardDestination(rhs.Checkerboard(), shift, dimension);
  Cshift_local(ret, rhs, dimension, shift);
  return ret;
}
NAMESPACE_END(Grid);

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
  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(&Grid);

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > coor(&Grid);
  for(int mu=0;mu<4;mu++) {
    LatticeComplex tmp(&Grid);
    LatticeCoordinate(tmp, mu);
    pokeLorentz(coor, tmp, mu);
  }

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > half_e(rbGrid);
  pickCheckerboard(Even, half_e, coor);
  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > half_o(rbGrid);
  pickCheckerboard(Odd, half_o, coor);

  // half_o = Cshift(half_o, 0, 1);
  // half_o.checkerboard = Even;

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > shifted(rbGrid);
  int dimension = 3, shift= 1;
  shifted = my_Cshift_local_rb(half_o, dimension, shift);
  // Cshift_local(shifted, half_o, dimension, shift);
  // shifted.Checkerboard() = Odd;  // Must manually set checkerboard

  print_half_field(shifted);

  // std::cout << shifted << std::endl;

  Grid_finalize();
}


