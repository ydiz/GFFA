#include <cmath>

namespace Grid{
namespace QCD{

void GF_heatbath(const LatticeGaugeField &Umu, LatticeColourMatrix &g,
                int nsweeps, Real _betaMM, const std::string &table_path,
                LatticeGaugeField *dSGF2dU=NULL,
                LatticeGaugeField (* const ForceFunc)(const LatticeColourMatrix &, const std::vector<LatticeColourMatrix> &)=NULL)//, bool verbose=0)
{
  static bool initialized = false;
  static GridRedBlackCartesian rbGrid(Umu._grid);
  static GridParallelRNG  pRNG(&rbGrid);
  static GridSerialRNG    sRNG;
  if (!initialized) {
   initialized = true;
   pRNG.SeedFixedIntegers(std::vector<int>{1,2,3,4});
   sRNG.SeedFixedIntegers(std::vector<int>{5,6,7,8});
  }

  RealD coeff = _betaMM;
  // Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  std::vector<LatticeColourMatrix> g_oe(2, &rbGrid);
  for(int cb=0; cb<2; cb++) pickCheckerboard(cb, g_oe[cb], g);

  std::vector<std::vector<LatticeColourMatrix>> U_oe(4, std::vector<LatticeColourMatrix>(2, &rbGrid));
  std::vector<std::vector<LatticeColourMatrix>> UMinusShift_oe(4, std::vector<LatticeColourMatrix>(2, &rbGrid));
  {
    LatticeColourMatrix UMinusShift_tmp(Umu._grid);
    for(int mu=0; mu<4; ++mu) {
      UMinusShift_tmp = Cshift(U[mu], mu, -1);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, U_oe[mu][cb], U[mu]);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, UMinusShift_oe[mu][cb], UMinusShift_tmp);
    }
  }

  LatticeColourMatrix staple_half(&rbGrid);
  // FIXME: can be optimized; only calculate staple_half!!!
  // most straightforard way: pickCheckerboard(Cshift(g,mu,1), cb)
  for(int sweep=0; sweep<nsweeps; sweep++){
  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for(int cb=0;cb<2;cb++) {

      staple_half = zero;
      staple_half.checkerboard = cb;
      int cb_inverse = (cb==1) ? 0 : 1;
      //if we choose the sign of S_GF1 to be postive, sign of staple should be negative
      for(int mu=0; mu<Nd; ++mu) {
        // Cshift will adjuct checkerboard automatically
        staple_half = staple_half + U_oe[mu][cb] * adj(Cshift(g_oe[cb_inverse],mu,1))
                        + adj( Cshift(g_oe[cb_inverse], mu, -1) * UMinusShift_oe[mu][cb]);
      }

      for(int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        GF_SubGroupHeatBath(sRNG, pRNG, coeff, g_oe[cb], staple_half, subgroup, cb, table_path);
      }

    }
    setCheckerboard(g, g_oe[Odd]);
    setCheckerboard(g, g_oe[Even]);  // FIXME: change ForceFunc so that we can operate only on g_oe.

    if(sweep%10==9) ProjectOnGroup(g); // project on group every 10 sweeps

    if(dSGF2dU!=NULL) (*dSGF2dU) += ForceFunc(g, U); // FIXME: is it better to have an inteval between added samples

  }



} // end of whole function




}}
