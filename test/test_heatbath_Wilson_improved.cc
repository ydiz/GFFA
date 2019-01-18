#include <Grid/Grid.h>
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
  // readField(Umu, "./U_wilson_equilbrium_8888");
  // readField(Umu, "./a0.2M3e0.2_ckpoint_lat");
  // Umu=1.0; // Cold start
  GridParallelRNG  pRNG_full(grid); pRNG_full.SeedFixedIntegers(std::vector<int>{1,2,3,4});
  SU3::HotConfiguration(pRNG_full, Umu);

  // RNG set up for test
  std::vector<int> pseeds({5,6,7,8,9}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(rbGrid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // SU3 colour operatoions
  LatticeColourMatrix link(grid);
  LatticeColourMatrix staple(grid);

  // Apply heatbath to the link
  // RealD beta=5.6;
  RealD beta = 6.0;
  RealD k_bar = 4.0;
  RealD coeff = beta / 3.0;
  std::string table_path = "/home/yz/GFFA/jupyter/numerical_integration/lookup_table_beta6.0";

  LatticeColourMatrix staple_half(rbGrid);
  LatticeColourMatrix U_half(rbGrid);

  for(int sweep=0;sweep<40;sweep++){

    RealD plaq = ColourWilsonLoops::avgPlaquette(Umu);

    std::cout<<GridLogMessage<<"sweep "<<sweep<<" PLAQUETTE "<<plaq<<std::endl;

    for( int cb=0;cb<2;cb++ ) {

      for(int mu=0;mu<Nd;mu++){

      	ColourWilsonLoops::Staple(staple,Umu,mu);
      	link = PeekIndex<LorentzIndex>(Umu,mu);

        pickCheckerboard(cb, staple_half, staple);
        pickCheckerboard(cb, U_half, link);

      	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {

          GF_SubGroupHeatBath(sRNG, pRNG, coeff, U_half, staple_half, subgroup, cb, table_path);

      	}

        setCheckerboard(link, U_half);

      	PokeIndex<LorentzIndex>(Umu,link,mu);
      	//reunitarise link;
      	// ProjectOnGroup(Umu);
      }
    }

  }

  Grid_finalize();
}
