namespace Grid {
namespace QCD {

//add betaMM

template <class Gimpl>
class GFAction : public Action<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  explicit GFAction(RealD beta_, RealD betaMM_, int innerMC_N_, int hb_offset_, const std::string& table_path_)
  : beta(beta_), betaMM(betaMM_), innerMC_N(innerMC_N_), hb_offset(hb_offset_), table_path(table_path_){}

  virtual std::string action_name() {return "GFAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[GFAction] Beta: " << beta  <<  "\t BetaMM: " << betaMM
      << "\t innerMC_N: " <<innerMC_N << "\t hb_offset: " << hb_offset <<  std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  //delta S_GF2 is calculated in GF_HMC.h: evolve();
  virtual RealD S(const GaugeField &U) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    RealD Sw = Waction.S(U);
  	std::cout  << "Wilson S: " <<  std::setprecision(15) << Sw << std::endl;

    RealD SGF1 = - (1./3.) * betaMM * Omega_no_g(U);
  	std::cout  << "SGF1: " <<  std::setprecision(15) << SGF1 << std::endl;

    return Sw + SGF1;
  }

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    GaugeField dSwdU(U.Grid());
    Waction.deriv(U, dSwdU);

  	RealD factor = 0.5 * (1./3.) * betaMM;

    GaugeField dSGF1dU(U.Grid());
    dSGF1dU = factor * Ta(U);

    // dSdU = dSwdU + dSGF1dU;

    GaugeField dSGF2dU(U.Grid());
    dSGF2dU = Zero();
    static LatticeColourMatrix g(U.Grid());
    static bool g_initialized = false;
    if(! g_initialized) {
      g = 1.0;
      g_initialized = true;
    }

    GF_heatbath(U, g, hb_offset, betaMM, table_path); //hb_nsweeps before calculate equilibrium value
    GF_heatbath(U, g, innerMC_N, betaMM, table_path, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU

    dSGF2dU = factor *  (1.0 / double(innerMC_N)) * dSGF2dU;

    dSdU = dSwdU + dSGF1dU - dSGF2dU;

  }
private:
  RealD beta;
  RealD betaMM;
  int innerMC_N;
  int hb_offset;
  std::string table_path;
};


template <class Gimpl>
class GF_DBW2Action : public Action<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  explicit GF_DBW2Action(RealD beta_, RealD betaMM_, int innerMC_N_, int hb_offset_, const std::string& table_path_)
  : beta(beta_), betaMM(betaMM_), innerMC_N(innerMC_N_), hb_offset(hb_offset_), table_path(table_path_){}

  virtual std::string action_name() {return "GF_DBW2Action";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[GF_DBW2Action] Beta: " << beta  <<  "\t BetaMM: " << betaMM
      << "\t innerMC_N: " <<innerMC_N << "\t hb_offset: " << hb_offset <<  std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  //delta S_GF2 is calculated in GF_HMC.h: evolve();
  virtual RealD S(const GaugeField &U) {
    DBW2GaugeAction<Gimpl> DBW2_action(beta);

    RealD Sw = DBW2_action.S(U);
  	std::cout  << "DBW2 S: " <<  std::setprecision(15) << Sw << std::endl;

    RealD SGF1 = - betaMM * Omega_no_g(U);
  	std::cout  << "SGF1: " <<  std::setprecision(15) << SGF1 << std::endl;

    return Sw + SGF1;
  }

  virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
    
    DBW2GaugeAction<Gimpl> DBW2_action(beta);
    GaugeField dSwdU(U.Grid());
    DBW2_action.deriv(U, dSwdU);

  	RealD factor = 0.5 * betaMM;

    GaugeField dSGF1dU(U.Grid());
    dSGF1dU = factor * Ta(U);

    GaugeField dSGF2dU(U.Grid());
    dSGF2dU = Zero();
    LatticeColourMatrix g(U.Grid());
    g = 1.0;
    GF_heatbath(U, g, hb_offset, betaMM, table_path); //hb_nsweeps before calculate equilibrium value

    GF_heatbath(U, g, innerMC_N, betaMM, table_path, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU

    dSGF2dU = factor *  (1.0 / double(innerMC_N)) * dSGF2dU;

    dSdU = dSwdU + dSGF1dU - dSGF2dU;

  }
private:
  RealD beta;
  RealD betaMM;
  int innerMC_N;
  int hb_offset;
  std::string table_path;
};

typedef GFAction<PeriodicGimplR>          GFActionR;
typedef GF_DBW2Action<PeriodicGimplR>     GF_DBW2ActionR;

}}
