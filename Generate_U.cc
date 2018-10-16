#include <Grid/Grid.h>
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

  GridRedBlackCartesian * rbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  LatticeGaugeField Umu(grid);
  Umu=1.0; // Cold start
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);

  // Apply heatbath to the link
  RealD beta=5.6;

  int subsets[2] = {Even, Odd};
  LatticeInteger one(rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(grid);
  LatticeColourMatrix link(grid);
  LatticeColourMatrix staple(grid);

  for(int sweep=0;sweep<2000;sweep++){
    RealD plaq = ColourWilsonLoops::avgPlaquette(Umu);
    std::cout<<GridLogMessage<<"sweep "<<sweep<<" PLAQUETTE "<<plaq<<std::endl;
    for( int cb=0;cb<2;cb++ ) {
      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);
      //      std::cout<<GridLogMessage<<mask<<std::endl;
      for(int mu=0;mu<Nd;mu++){
      	// Get Link and Staple term in action; must contain Beta and any other coeffs
      	ColourWilsonLoops::Staple(staple,Umu,mu);
      	link = PeekIndex<LorentzIndex>(Umu,mu);
      	for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++ ) {
      	  // update Even checkerboard
      	  SU3::SubGroupHeatBath(sRNG,pRNG,beta,link,staple,subgroup,20,mask);
      	}
      	PokeIndex<LorentzIndex>(Umu,link,mu);
      	//reunitarise link;
      	ProjectOnGroup(Umu);
      }
    }
  }

  LatticeColourMatrix g(Umu._grid);
  g = 1.0;
  double M = 1.0;
  double betaMM = beta * M * M;
  GF_heatbath(Umu, g, 10000, betaMM, 10);
  SU3::GaugeTransform(Umu, g);
  // std::cout << "after SDGF dOmegaSquare: " << dOmegaSquare2(g, U) << std::endl;
  std::string fileName="U8x8_M1.0_beta5.6";
  NerscIO::writeConfiguration(Umu, fileName, 0, 0);

  Grid_finalize();
}
