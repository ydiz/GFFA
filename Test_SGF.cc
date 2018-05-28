#include <Grid/Grid.h>
#include <iostream>
#include "GF_Util.h"
#include "GF_heatbath_Util.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  std::vector<int> latt({8,8,8,8});
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(latt,
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  // RNG set up for test
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  LatticeGaugeField U(grid);
  FieldMetaData header;
  std::string file("./U_equilibrium");
  NerscIO::readConfiguration(U, header,file);

  LatticeColourMatrix g(grid);
  g = 1.0;

  double betaMM = 50;
  int sweeps = 5000;
  bool SDGF = true;

  stringstream ss;
  string arg;

  if(GridCmdOptionExists(argv, argv+argc, string("--betaMM"))){
    arg = GridCmdOptionPayload(argv, argv+argc, string("--betaMM"));
    ss.clear();
    ss.str(arg);
    ss >> betaMM;
  }

  if(GridCmdOptionExists(argv, argv+argc, string("--sweeps"))){
    arg = GridCmdOptionPayload(argv, argv+argc, string("--sweeps"));
    ss.clear();
    ss.str(arg);
    ss >> sweeps;
  }

  if(GridCmdOptionExists(argv, argv+argc, string("--SDGF"))){
    arg = GridCmdOptionPayload(argv, argv+argc, string("--SDGF"));
    ss.clear();
    ss.str(arg);
    ss >> SDGF;
  }

  cout << "sweeps: " << sweeps << endl;
  cout << "betaMM: " << betaMM << endl;
  cout << "SDGF: " << SDGF << endl;

  Real alpha=0.1;
  // must be false. cannot be true
  // if(SDGF) FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U,alpha,200,1.0e-12, 1.0e-12,false);

  if(SDGF){
    Real alpha=0.1;
    // must be false. cannot be true
    std::cout << "before SDGF dOmegaSquare: " << dOmegaSquare2_no_g(U) << std::endl;
    std::cout << "before SDGF SGF1: " << - betaMM * Omega_no_g(U) << std::endl;

    LatticeColourMatrix g(U._grid);
    g = 1.0;
    GF_heatbath(U, g, 100, betaMM);
    SU3::GaugeTransform(U, g);
    //FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U,alpha,100,1.0e-12, 1.0e-12,false);
    std::cout << "after SDGF dOmegaSquare: " << dOmegaSquare2_no_g(U) << std::endl;
    std::cout << "after SDGF SGF1: " << - betaMM * Omega_no_g(U) << std::endl;
  }

  GF_heatbath(U, g, sweeps, betaMM, 1);




  Grid_finalize();
} // main
