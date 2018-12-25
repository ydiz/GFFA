#include <Grid/Grid.h>
// #include "../subgroup_hb_rbgrid.h"
#include "../subgroup_hb.h"

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

  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  ///////////////////////////////
  // Configuration of known size
  ///////////////////////////////
  GridParallelRNG  pRNG_full(grid); pRNG_full.SeedFixedIntegers(std::vector<int>{1,2,3,4});
  LatticeGaugeField Umu(grid);
  SU3::HotConfiguration(pRNG_full, Umu);

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  std::vector<LatticeColourMatrix> UMinusShift(4, Umu._grid);
  for(int mu=0; mu<Nd; ++mu) UMinusShift[mu] = Cshift(U[mu], mu, -1);


  // RNG set up for test
  std::vector<int> pseeds({5,6,7,8,9}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(rbGrid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // SU3 colour operatoions
  LatticeColourMatrix g(grid);
  g = 1.0;
  LatticeColourMatrix staple(grid);

  // Apply heatbath to the link
  RealD beta=5.6;

  LatticeColourMatrix staple_half(rbGrid);
  LatticeColourMatrix g_half(rbGrid);

  for(int sweep=0;sweep<40;sweep++){
    // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;

    for( int cb=0;cb<2;cb++ ) {
      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * UMinusShift[mu]);
      }

        pickCheckerboard(cb, staple_half, staple);
        pickCheckerboard(cb, g_half, g);

      	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {
          GF_SubGroupHeatBath(sRNG,pRNG,beta,g_half,staple_half,subgroup,1, cb);
      	}

        setCheckerboard(g, g_half);

    }

    if(sweep == 2) {
      cout << g << endl;
      assert(0);
    }

  }

  Grid_finalize();
}
