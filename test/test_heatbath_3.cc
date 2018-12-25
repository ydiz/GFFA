#include <Grid/Grid.h>
// #include "../subgroup_hb_rbgrid.h"
#include "../subgroup_hb.h"
#include "../GF_Util.h"

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

  // RNG set up for test
  std::vector<int> pseeds({5,6,7,8,9}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(rbGrid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // SU3 colour operatoions
  LatticeColourMatrix g(grid);
  g = 1.0;
  // LatticeColourMatrix staple(grid);

  // Apply heatbath to the link
  RealD beta=5.6;

  // LatticeColourMatrix g_half(rbGrid);
  std::vector<LatticeColourMatrix> g_oe(2, rbGrid);
  for(int cb=0; cb<2; cb++) pickCheckerboard(cb, g_oe[cb], g);

  std::vector<std::vector<LatticeColourMatrix>> U_oe(4, std::vector<LatticeColourMatrix>(2, rbGrid));
  std::vector<std::vector<LatticeColourMatrix>> UMinusShift_oe(4, std::vector<LatticeColourMatrix>(2, rbGrid));
  {
    LatticeColourMatrix U_tmp(Umu._grid);
    LatticeColourMatrix UMinusShift_tmp(Umu._grid);
    for(int mu=0; mu<4; ++mu) {
      U_tmp = PeekIndex<LorentzIndex>(Umu, mu);
      UMinusShift_tmp = Cshift(U_tmp, mu, -1);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, U_oe[mu][cb], U_tmp);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, UMinusShift_oe[mu][cb], UMinusShift_tmp);
    }
  }

  LatticeColourMatrix staple_half(rbGrid);
  for(int sweep=0;sweep<40;sweep++){
    // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for( int cb=0;cb<2;cb++ ) {
      staple_half = zero;
      staple_half.checkerboard = cb;
      int cb_inverse = (cb==1) ? 0 : 1;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        // staple += U[mu] * adj(Cshift(g,mu,1)) + adj( Cshift(g, mu, -1) * UMinusShift[mu]);
      staple_half = staple_half + U_oe[mu][cb] * adj(Cshift(g_oe[cb_inverse],mu,1))
                      + adj( Cshift(g_oe[cb_inverse], mu, -1) * UMinusShift_oe[mu][cb]);
      }

    	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {
        GF_SubGroupHeatBath(sRNG, pRNG, beta, g_oe[cb], staple_half, subgroup, 1, cb);
    	}
    }

    if(sweep == 2) {
      setCheckerboard(g, g_oe[Odd]);
      setCheckerboard(g, g_oe[Even]);
      cout << g << endl;
      assert(0);
    }

  }

  Grid_finalize();
}
