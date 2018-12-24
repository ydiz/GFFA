#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  // std::vector<int> latt({8,8,8,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > coor(grid);
  LatticeComplex tmp(grid);

  for(int mu=0;mu<4;mu++){
    LatticeCoordinate(tmp, mu);
    pokeLorentz(coor, tmp, mu);
  }

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > half_e(rbGrid);
  pickCheckerboard(Even, half_e, coor);

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > half_o(rbGrid);
  pickCheckerboard(Odd, half_o, coor);

  half_o = Cshift(half_o, 0, 1);
  half_o.checkerboard = Even;

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > tmp_half = half_o - half_e;

  Lattice< iVector<iScalar<iScalar<vComplex>>, Nd> > output(grid);
  output = zero;


  // auto tmp_half = Cshift(half, 0, 1);
  setCheckerboard(output, tmp_half);

  cout << output << endl;


  Grid_finalize();
}
