#include <Grid/Grid.h>
#include <iostream>

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

  // RNG set up for test
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // LatticeGaugeField U(grid);
  // FieldMetaData header;
  // std::string file("./U_equilibrium");
  // NerscIO::readConfiguration(U, header,file);

  LatticeGaugeField P(grid);
  PeriodicGimplR::generate_momenta(P, pRNG);
  std::string file("P_4by4");
  NerscIO::writeConfiguration(P,file,0,0);

  // FieldMetaData header;
  // std::string file("./P_4by4");
  // NerscIO::readConfiguration(P, header, file);

  cout << P;

}
