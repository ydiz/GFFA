#include <Grid/Grid.h>
#include "../GF_Util.h"
#include "../GF_init_k.h"
#include "../GF_generate_P.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);

  LatticeGaugeField P(grid);
  // PeriodicGimplR::generate_momenta(P, pRNG);
  double M = 1.0, epsilon = 0.2;
  if(argc >=2) M = std::stod(std::string(argv[1]));
  if(argc >=3) epsilon = std::stod(std::string(argv[2]));

  const Momenta_k KK(P._grid, M, epsilon, 0);
  GF_generate_P(P, pRNG, KK);

  // calculate average and standard deviation
  std::vector<ColourMatrix> ta(8);
  for(int i=0; i<8; ++i) SU3::generator(i, ta[i]);

  LatticeGaugeField tt(P._grid);
  Lattice<iVector<iScalar<iScalar<vReal > > ,4>> tmp(P._grid);

  for(int a=0; a<8; ++a)
  {
    tmp = toReal(timesI(2.0 * trace(P * ta[a])));
    auto summation = sum(tmp);
    double avg = ((summation(0) + summation(1) + summation(2) + summation(3))()()) / P._grid -> gSites() / 4.0;
    tmp = tmp * tmp;
    auto tmp_sum = sum(tmp) ; // tmp is vector<scalar<scalar<Real>>>>
    double variance = ((tmp_sum(0) + tmp_sum(1) + tmp_sum(2) + tmp_sum(3))()()) / P._grid -> gSites() / 4.0 ;
    double std_deviation = std::sqrt(variance);
    std::cout << "a = " << a << ": " << "average: "<< avg << " variance: "<< variance << " std: " << std_deviation << std::endl;
  }

  Grid_finalize();
}
