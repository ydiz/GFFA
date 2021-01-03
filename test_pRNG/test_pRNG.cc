#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  GridParallelRNG pRNG(grid);  // not used

  pRNG.SeedFixedIntegers(std::vector<int>({10,20,30,40}));
  
  LatticeComplex lat(grid);
  random(pRNG, lat);
  std::cout << lat << std::endl;

  std::cout << "Finished!"  << std::endl;
  Grid_finalize();

} // main
