#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({4,4,4,4});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt,
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  int seed = std::chrono::system_clock::now().time_since_epoch().count();

  GridParallelRNG pRNG(grid);
  pRNG.SeedFixedIntegers(std::vector<int>{seed});
  LatticeGaugeField U(grid);
  gaussian(pRNG, U);

  cout << U << endl;

     cout << seed << endl; //1646588203

}
