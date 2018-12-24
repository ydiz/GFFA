#include <Grid/Grid.h>
#include "../subgroup_hb_rbgrid.h"
// #include "../subgroup_hb.h"

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
  LatticeGaugeField Umu(grid);
  Umu=1.0; // Cold start

  // RNG set up for test
  std::vector<int> pseeds({5,6,7,8,9}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(rbGrid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // SU3 colour operatoions
  LatticeColourMatrix link(grid);
  LatticeColourMatrix staple(grid);

  // Apply heatbath to the link
  RealD beta=5.6;

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
      	  // update Even checkerboard
          GF_SubGroupHeatBath(sRNG,pRNG,beta,U_half,staple_half,subgroup,1, cb);
          // GF_SubGroupHeatBath(sRNG,pRNG,beta,U_half,staple_half,subgroup,20, cb);
          // assert(0);
      	}

        // std::cout << "U_half checkerboard: " << U_half.checkerboard << "\n";

        setCheckerboard(link, U_half);


      	PokeIndex<LorentzIndex>(Umu,link,mu);
      	//reunitarise link;
      	// ProjectOnGroup(Umu);
      }



    }



  }

  // std::cout << Umu << std::endl;
  // assert(0);






  Grid_finalize();
}
