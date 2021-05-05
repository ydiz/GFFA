#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  GridSerialRNG sRNG;

  sRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  for(int i=0; i<10; ++i) {
    iMatrix<Complex, 3> m;
    random(sRNG, m);
    std::cout << "before ProjectOnGroup" << std::endl;
    std::cout << m << std::endl;
    m = ProjectOnGroup(m);
    std::cout << "after ProjectOnGroup" << std::endl;
    std::cout << m << std::endl;
    std::cout << "after ProjectOnGroup adj(m) * m" << std::endl;
    std::cout << adj(m) * m << std::endl;
    std::cout << "===========================================" << std::endl;
  }
  
  // LatticeComplex lat(grid);
  // random(pRNG, lat);
  // std::cout << lat << std::endl;

  std::cout << "Finished!"  << std::endl;
  Grid_finalize();

} // main
