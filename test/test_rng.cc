#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({8,8,8,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  std::vector<int> pseeds({5,6,7,8,9}); // once I caught a fish alive
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridParallelRNG  pRNG_rb(rbGrid); pRNG_rb.SeedFixedIntegers(std::vector<int>{1,2,3,4,13});


  LatticeReal cos_theta(grid);
  LatticeReal tmp_half(rbGrid);
  tmp_half.checkerboard = 0;

  cos_theta = 0.0;
  random(pRNG, cos_theta);
  std::cout << cos_theta << std::endl;
  std::cout << sum(cos_theta) / cos_theta._grid->gSites() << std::endl;

  // cos_theta = 0.0;
  // random(pRNG_rb, tmp_half);
  // setCheckerboard(cos_theta, tmp_half);
  // // std::cout << cos_theta << std::endl;
  // std::cout << "average: "<< sum(cos_theta) / tmp_half._grid->gSites() << std::endl;
  // std::cout << "grid: "<< rbGrid->iSites() << " " << rbGrid->oSites() << " lSites: "<< rbGrid->lSites() << std::endl;

}
