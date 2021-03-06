#include <cmath>

namespace Grid{

void GF_heatbath(const LatticeGaugeField &Umu, LatticeColourMatrix &g,
                int nsweeps, Real _betaMM, const std::string &table_path,
                GridParallelRNG &pRNG,
                LatticeGaugeField *dSGF2dU=NULL,
                LatticeGaugeField (* const ForceFunc)(const LatticeColourMatrix &, const std::vector<LatticeColourMatrix> &)=NULL)//, bool verbose=0)
{

  GridRedBlackCartesian rbGrid(Umu.Grid());

  RealD coeff = _betaMM * (1./3.);
  // Real betaMM = _betaMM * 3; //In heatbath routine, the coefficient is beta/Nc

  std::vector<LatticeColourMatrix> U(4, Umu.Grid());
  for(int mu=0; mu<Nd; mu++) U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  std::vector<LatticeColourMatrix> g_oe(2, &rbGrid);
  for(int cb=0; cb<2; cb++) pickCheckerboard(cb, g_oe[cb], g);

  std::vector<std::vector<LatticeColourMatrix>> U_oe(4, std::vector<LatticeColourMatrix>(2, &rbGrid));
  std::vector<std::vector<LatticeColourMatrix>> UMinusShift_oe(4, std::vector<LatticeColourMatrix>(2, &rbGrid));
  {
    LatticeColourMatrix UMinusShift_tmp(Umu.Grid());
    for(int mu=0; mu<4; ++mu) {
      UMinusShift_tmp = Cshift(U[mu], mu, -1);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, U_oe[mu][cb], U[mu]);
      for(int cb=0; cb<2; cb++) pickCheckerboard(cb, UMinusShift_oe[mu][cb], UMinusShift_tmp);
    }
  }

  std::vector<std::vector<LatticeColourMatrix>> dSGF2dU_oe(4, std::vector<LatticeColourMatrix>(2, &rbGrid));
  if(dSGF2dU!=NULL) {
    for(int cb=0;cb<2;cb++) {
      for(int mu=0;mu<4;mu++) {
        dSGF2dU_oe[mu][cb] = Zero();
        dSGF2dU_oe[mu][cb].Checkerboard() = cb;
      }
    }
  }

  for(int sweep=0; sweep<nsweeps; sweep++){
  // std::cout<<GridLogMessage<<"sweep "<<sweep<<" S: "<<GF_S(g, Umu)<<std::endl;
    for(int cb=0;cb<2;cb++) {

      LatticeColourMatrix staple_half(&rbGrid);
      staple_half = Zero();
      staple_half.Checkerboard() = cb;
      int cb_inverse = (cb==1) ? 0 : 1;
      
      for(int mu=0; mu<Nd; ++mu) {
        // Cshift will adjust Checkerboard() automatically
        // if we choose the sign of S_GF1 to be postive, sign of staple should be negative
        staple_half += U_oe[mu][cb] * adj(Cshift(g_oe[cb_inverse],mu,1))
                        + adj( Cshift(g_oe[cb_inverse], mu, -1) * UMinusShift_oe[mu][cb]);
      }

      for(int subgroup=0;subgroup<SU3::su2subgroups();subgroup++) {
        GF_SubGroupHeatBath(pRNG, coeff, g_oe[cb], staple_half, subgroup, cb, table_path);
      }


    }



    // I found that project on group has almost no effect. at least on a small lattice
    // if(sweep%10==9) {ProjectOnGroup(g_oe[Odd]); ProjectOnGroup(g_oe[Even]);}// project on group every 10 sweeps

    if(dSGF2dU!=NULL) {
      for(int cb=0;cb<2;cb++) {
        int cb_inverse = (cb==1) ? 0 : 1;
        for(int mu=0; mu<Nd; ++mu) {
          dSGF2dU_oe[mu][cb] += U_oe[mu][cb] * adj(Cshift(g_oe[cb_inverse], mu, 1)) * g_oe[cb]; //FIXME: if calculate force every sweep, I can save Csfhit g
        }
  	  }
  	}

    // if(dSGF2dU!=NULL) (*dSGF2dU) += ForceFunc(g, U); // FIXME: is it better to have an inteval between added samples
  } // end of sweeps loop

  setCheckerboard(g, g_oe[Odd]);
  setCheckerboard(g, g_oe[Even]);

  if(dSGF2dU!=NULL) {
    LatticeColourMatrix tmp(Umu.Grid());
    for(int mu=0;mu<4;mu++) {
        setCheckerboard(tmp, dSGF2dU_oe[mu][Odd]);
        setCheckerboard(tmp, dSGF2dU_oe[mu][Even]);
        PokeIndex<LorentzIndex>(*dSGF2dU, tmp, mu);
    }
    *dSGF2dU = Ta(*dSGF2dU);
  }
  

  if(dSGF2dU!=NULL) {
    // // print ||adj(g) * g - I||^2
    // LatticeColourMatrix one(Umu.Grid());
    // one = 1.0;
    // LatticeColourMatrix tmp = adj(g) * g - one;
    // std::cout << "||adj(g) * g - I||^2: " << norm2(tmp) << std::endl;

    // Re-unitarize Umu
    g = ProjectOnGroup(g);
    // tmp = adj(g) * g - one;
    // std::cout << "||adj(g) * g - I||^2" << norm2(tmp) << std::endl;
  }



} // end of whole function




}
