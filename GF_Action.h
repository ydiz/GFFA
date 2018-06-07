namespace Grid {
namespace QCD {

//add betaMM

template <class Gimpl>
class GFAction : public Action<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  /////////////////////////// constructors
  explicit GFAction(RealD beta_, RealD betaMM_, int innerMC_N_, int hb_offset_, int hb_nsweeps_, int hb_multi_hit_): beta(beta_), betaMM(betaMM_), innerMC_N(innerMC_N_), hb_offset(hb_offset_), hb_nsweeps(hb_nsweeps_), hb_multi_hit(hb_multi_hit_){}

  virtual std::string action_name() {return "GFAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[GFAction] Beta: " << beta  <<  "\t BetaMM: " << betaMM
      << "\t innerMC_N: " <<innerMC_N << "\t hb_offset: " << hb_offset << "\t hb_nsweeps: " <<hb_nsweeps << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  //delta S_GF2 is calculated in GF_HMC.h: evolve();
  virtual RealD S(const GaugeField &U) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    RealD Sw = Waction.S(U);
  	std::cout << "Wilson S: " << Sw << std::endl;

    RealD SGF1 = - betaMM * Omega_no_g(U);
  	std::cout << "SGF1: " << SGF1 << std::endl;

    return Sw + SGF1;
  }

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    GaugeField dSwdU(U._grid);
    Waction.deriv(U, dSwdU);

  	RealD factor = 0.5 * betaMM;

    GaugeField dSGF1dU(U._grid);
    dSGF1dU = factor * Ta(U);

    GaugeField dSGF2dU(U._grid);
    dSGF2dU = zero;
    LatticeColourMatrix g(U._grid);
    g = 1.0;
    GF_heatbath(U, g, hb_offset, betaMM, hb_multi_hit); //hb_nsweeps before calculate equilibrium value

    //???maybe error here is large;(if dSGF2dU is large)
    for(int i=0; i<innerMC_N; ++i)
    {
      dSGF2dU += dOmegadU_g(g, U);
      GF_heatbath(U, g, hb_nsweeps, betaMM, hb_multi_hit);
    }

    dSGF2dU = factor *  (1.0 / double(innerMC_N)) * dSGF2dU;

    dSdU = dSwdU + dSGF1dU - dSGF2dU;

  }
private:
  RealD beta;
  RealD betaMM;
  int innerMC_N;
  int hb_offset;
  int hb_nsweeps;
  int hb_multi_hit;
};


typedef GFAction<PeriodicGimplR>          GFActionR;

}}
