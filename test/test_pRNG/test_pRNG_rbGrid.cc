#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  GridParallelRNG pRNG(rbGrid);
  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));
  
  LatticeComplex lat(rbGrid);
  random(pRNG, lat);

  autoView(lat_v, lat, AcceleratorRead);
  for(int ss=0; ss<lat.Grid()->oSites(); ++ss) {
    std::cout << ss << " ";
    std::cout << lat_v[ss] << std::endl;
  }

  std::cout << "Finished!"  << std::endl;
  Grid_finalize();

} // main
