#include <cmath>

namespace Grid{
namespace QCD{

//betaMM; use pRNG, sRNG outside
void GF_heatbath(const LatticeGaugeField &Umu, LatticeColourMatrix &g,
                int nsweeps, Real _betaMM, int multi_hit,
                LatticeGaugeField *dSGF2dU=NULL,
                LatticeGaugeField (*ForceFunc)(const LatticeColourMatrix &, const std::vector<LatticeColourMatrix> &)=NULL)//, bool verbose=0)
{
  GridRedBlackCartesian rbGrid(Umu._grid);

  int seed = 1;
  // int seed = std::chrono::system_clock::now().time_since_epoch().count();
  //FIXME: seeds are the same for each step
  GridParallelRNG  pRNG(Umu._grid); pRNG.SeedFixedIntegers(std::vector<int>{seed});
  GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(std::vector<int>{seed});

  Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  int subsets[2] = {Even, Odd};
  LatticeInteger one(&rbGrid);  one = 1; // fill with ones
  LatticeInteger mask(Umu._grid);

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix staple(Umu._grid);

  for(int sweep=0; sweep<nsweeps; sweep++){
  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for(int cb=0;cb<2;cb++) {
      one.checkerboard=subsets[cb];
      mask= zero;
      setCheckerboard(mask,one);

      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * Cshift(U[mu], mu, -1) ); // FIXME: calculate and save Cshift(U[mu], mu, -1)
      }
      for(int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        SU3::SubGroupHeatBath(sRNG,pRNG,betaMM,g,staple,subgroup,multi_hit,mask);
      }
      //reunitarise
      ProjectOnGroup(g);
    }

  if(dSGF2dU!=NULL) (*dSGF2dU) += ForceFunc(g, U);

  }

}




}}
