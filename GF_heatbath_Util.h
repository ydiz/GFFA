#include <cmath>

namespace Grid{
namespace QCD{

// FIXME: "one" and "mask" can be defined outside and generated only once.
void GF_heatbath(const LatticeGaugeField &Umu, LatticeColourMatrix &g,
                int nsweeps, Real _betaMM, int multi_hit,
                LatticeGaugeField *dSGF2dU=NULL,
                LatticeGaugeField (* const ForceFunc)(const LatticeColourMatrix &, const std::vector<LatticeColourMatrix> &)=NULL)//, bool verbose=0)
{
  static bool initialized = false;
  static GridParallelRNG  pRNG(Umu._grid);
  static GridSerialRNG    sRNG;
  if (!initialized) {
   initialized = true;
   // std::cout << "*************************" << "\n";
   pRNG.SeedFixedIntegers(std::vector<int>{1,2,3,4});
   sRNG.SeedFixedIntegers(std::vector<int>{5,6,7,8});
   // do the initialization part
}
  // int seed = 1;
  // int seed = std::chrono::system_clock::now().time_since_epoch().count();
  // GridParallelRNG  pRNG(Umu._grid); pRNG.SeedFixedIntegers(std::vector<int>{seed});
  // GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(std::vector<int>{seed});

  Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  GridRedBlackCartesian rbGrid(Umu._grid);
  int subsets[2] = {Even, Odd};
  // LatticeInteger one(&rbGrid);  one = 1; // fill with ones
  // LatticeInteger mask(Umu._grid);
  // added by zyd
  LatticeInteger one(&rbGrid); one = 1;
  std::vector<LatticeInteger> mask_EvenOdd(2, Umu._grid);
  for(int cb=0;cb<2;cb++) {
    one.checkerboard = subsets[cb];
    mask_EvenOdd[cb] = zero;
    setCheckerboard(mask_EvenOdd[cb], one);
  }

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix staple(Umu._grid);
  // added by zyd
  std::vector<LatticeColourMatrix> UMinusShift(4, Umu._grid);
  for(int mu=0; mu<Nd; ++mu) UMinusShift[mu] = Cshift(U[mu], mu, -1);

  for(int sweep=0; sweep<nsweeps; sweep++){
  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for(int cb=0;cb<2;cb++) {
      // one.checkerboard=subsets[cb];
      // mask= zero;
      // setCheckerboard(mask,one);

      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        // staple += U[mu] * adj(Cshift(g,mu,1))
        //           + adj( Cshift(g, mu, -1) * Cshift(U[mu], mu, -1) );
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * UMinusShift[mu]);
      }
      for(int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        // SU3::SubGroupHeatBath(sRNG,pRNG,betaMM,g,staple,subgroup,multi_hit,mask_EvenOdd[cb]);
        GF_SubGroupHeatBath(sRNG, pRNG, betaMM, g, staple, subgroup, multi_hit, mask_EvenOdd[cb]);
      }
      //reunitarise
	  ProjectOnGroup(g);
    }

    // if(sweep%10==9) ProjectOnGroup(g); // project on group every 10 sweeps

    if(dSGF2dU!=NULL) (*dSGF2dU) += ForceFunc(g, U);

  }

}




}}
