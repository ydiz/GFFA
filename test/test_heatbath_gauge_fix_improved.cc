#include <Grid/Grid.h>
// #include "../subgroup_hb_rbgrid.h"
#include "../Integral_table.h"
#include "../subgroup_hb.h"
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
  LatticeGaugeField Umu(grid);
  // readField(Umu, "./a0.2M3e0.2_ckpoint_lat");
  readField(Umu, "./U_wilson_equilbrium_8888");

  // GridParallelRNG  pRNG_full(grid); pRNG_full.SeedFixedIntegers(std::vector<int>{1,2,3,4});
  // LatticeGaugeField Umu(grid);
  // SU3::HotConfiguration(pRNG_full, Umu);

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
  // RealD k_bar = 6.457;
  // RealD beta = 6.0;
  RealD beta = 0.7796;
  RealD M = 4.0;
  RealD coeff = beta * M * M;
  std::string table_path = "/home/yz/GFFA/jupyter/numerical_integration/lookup_table_M4_beta0.7796";

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
    setCheckerboard(g, g_oe[Odd]);
    setCheckerboard(g, g_oe[Even]);
    std::cout<<GridLogMessage<<"sweep: "<< sweep <<" Omega_g: "<< Omega_g(g, Umu) << std::endl;
    std::cout<<GridLogMessage<<"sweep: "<<sweep<<" dOmegaSquare2: "<< dOmegaSquare2(g, Umu) << std::endl;
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
        GF_SubGroupHeatBath(sRNG, pRNG, coeff, g_oe[cb], staple_half, subgroup, 1, cb, table_path);
    	}
    }

  }

  Grid_finalize();
}
