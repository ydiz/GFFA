#include <Grid/Grid.h>
#include "../GF_Util.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
							GridDefaultSimd(Nd,vComplex::Nsimd()),
							GridDefaultMpi());

  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  ///////////////////////////////
  // Configuration of known size
  ///////////////////////////////
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  LatticeGaugeField Umu(grid);
  readField(Umu, "./U_wilson_equilbrium_8888");
  // LatticeGaugeField Umu(grid);
  // SU3::HotConfiguration(pRNG, Umu);
  // std::cout << Umu << std::endl; assert(0);

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  std::vector<LatticeColourMatrix> UMinusShift(4, Umu._grid);
  for(int mu=0; mu<Nd; ++mu) UMinusShift[mu] = Cshift(U[mu], mu, -1);


  // SU3 colour operatoions
  LatticeColourMatrix g(grid);
  g = 1.0;
  LatticeColourMatrix staple(grid);

  // Apply heatbath to the link
  RealD beta=6.0;
  RealD M = 3.0;
  RealD coeff = 3.0 * beta  * M * M;

  int subsets[2] = {Even, Odd};
  LatticeInteger one(rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(grid);

  for(int sweep=0;sweep<40;sweep++){
    std::cout<<GridLogMessage<<"sweep: "<<sweep<<" Omega_g: "<< Omega_g(g, Umu) << std::endl;
    std::cout<<GridLogMessage<<"sweep: "<<sweep<<" dOmegaSquare2: "<< dOmegaSquare2(g, Umu) << std::endl;

    for( int cb=0;cb<2;cb++ ) {
      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);

      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * UMinusShift[mu]);
      }

    	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {
        SU3::SubGroupHeatBath(sRNG,pRNG,coeff,g,staple,subgroup,1, mask);
    	}

    }

    // if(sweep == 2) {
    //   cout << g << endl;
    //   assert(0);
    // }

  }

  Grid_finalize();
}
