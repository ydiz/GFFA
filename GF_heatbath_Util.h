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
  static GridRedBlackCartesian rbGrid(Umu._grid);
  static GridParallelRNG  pRNG(&rbGrid);
  static GridSerialRNG    sRNG;
  if (!initialized) {
   initialized = true;
   // std::cout << "*************************" << "\n";
   pRNG.SeedFixedIntegers(std::vector<int>{1,2,3,4});
   sRNG.SeedFixedIntegers(std::vector<int>{5,6,7,8});
   // do the initialization part
}

  Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix staple(Umu._grid);
  // added by zyd
  std::vector<LatticeColourMatrix> UMinusShift(4, Umu._grid);
  for(int mu=0; mu<Nd; ++mu) UMinusShift[mu] = Cshift(U[mu], mu, -1);

  LatticeColourMatrix staple_half(&rbGrid);
  LatticeColourMatrix g_half(&rbGrid);

  // FIXME: can be optimized; only calculate staple_half!!!
  // most straightforard way: pickCheckerboard(Cshift(g,mu,1), cb)
  for(int sweep=0; sweep<nsweeps; sweep++){
  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for(int cb=0;cb<2;cb++) {
      staple = zero;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        staple += U[mu] * adj(Cshift(g,mu,1))
                  + adj( Cshift(g, mu, -1) * UMinusShift[mu]);
      }

      pickCheckerboard(cb, staple_half, staple);
      pickCheckerboard(cb, g_half, g);

      for(int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        // SU3::SubGroupHeatBath(sRNG,pRNG,betaMM,g,staple,subgroup,multi_hit,mask_EvenOdd[cb]);
        // GF_SubGroupHeatBath(sRNG, pRNG, betaMM, g_half, staple_half, subgroup, multi_hit, mask_EvenOdd[cb]);
        GF_SubGroupHeatBath(sRNG, pRNG, betaMM, g_half, staple_half, subgroup, multi_hit, cb);
      }

      setCheckerboard(g, g_half);
      //reunitarise
  	  // ProjectOnGroup(g);
    }

    if(sweep%10==9) ProjectOnGroup(g); // project on group every 10 sweeps

    if(dSGF2dU!=NULL) (*dSGF2dU) += ForceFunc(g, U);

  }

}




}}
