#include <cmath>

namespace Grid{
namespace QCD{

//betaMM; use pRNG, sRNG outside
void GF_heatbath(const LatticeGaugeField &Umu, LatticeColourMatrix &g, int nsweeps, Real _betaMM, bool verbose=0)
{
  GridRedBlackCartesian rbGrid(Umu._grid);

  //FIXME: seeds are the same for each step
  std::vector<int> pseeds({1,2,3,4,5}); // once I caught a fish alive
  std::vector<int> sseeds({6,7,8,9,10});// then i let it go again
  GridParallelRNG  pRNG(Umu._grid); pRNG.SeedFixedIntegers(pseeds);
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(sseeds);


  LatticeColourMatrix staple(Umu._grid);

  Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  int subsets[2] = { Even, Odd};
  LatticeInteger one(&rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(Umu._grid);

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  //dU_for_verbose and ta are for verbose; can be deleted
  //LatticeGaugeField dU_for_verbose(Umu._grid);
  //std::vector<ColourMatrix> ta(8);
  //for(int i=0; i<8; ++i) SU3::generator(i, ta[i]);

  for(int sweep=0;sweep<nsweeps;sweep++){

  //  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;

  //  if(verbose)
  //  {
  //    std::cout << GridLogMessage
  //    << "sweeps: " << sweep << std::endl;

  //    std::cout << GridLogMessage
  //    << "dOmegaSquare: " << dOmegaSquare2(g, Umu) << std::endl;

  //    dU_for_verbose = dSGF1dU_g(g, Umu);

  //    std::cout << GridLogMessage
  //    << "dU[1,2,3,4], 0: " << peek_xa(dU_for_verbose, std::vector<int>{1,2,3,4}, 0, ta) << std::endl;


  //    std::cout << GridLogMessage
  //    << "dU[2,2,2,2], 2" << peek_xa(dU_for_verbose, std::vector<int>{2,2,2,2}, 2, ta) << std::endl;

  //  }

    for( int cb=0;cb<2;cb++ ) {

      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);

      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu)
      {
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * Cshift(U[mu], mu, -1) );
      }

      for( int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        SU3::SubGroupHeatBath(sRNG,pRNG,betaMM,g,staple,subgroup,20,mask);
      }
      //reunitarise
      ProjectOnGroup(g);
    }

  }

}




}}
