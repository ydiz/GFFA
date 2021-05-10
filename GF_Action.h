namespace Grid {



template <class GaugeField >
class MyAction 
{

public:
  bool is_smeared = false;
  // Heatbath?
  virtual void refresh(const GaugeField& U, GridParallelRNG& pRNG) = 0; // refresh pseudofermions
  virtual RealD S(const GaugeField& U) = 0;                             // evaluate the action
  virtual void deriv(const GaugeField& U, GaugeField& dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) = 0;        // evaluate the action derivative
  virtual std::string action_name()    = 0;                             // return the action name
  virtual std::string LogParameters()  = 0;                             // prints action parameters
  virtual ~MyAction(){}
};



template <class Field, class Repr = NoHirep >
struct MyActionLevel {
public:
  unsigned int multiplier;

  // Fundamental repr actions separated because of the smearing
  typedef MyAction<Field>* ActPtr;

  // construct a tuple of vectors of the actions for the corresponding higher
  // representation fields
  typedef typename AccessTypes<MyAction, Repr>::VectorCollection action_collection;
  typedef typename  AccessTypes<MyAction, Repr>::FieldTypeCollection action_hirep_types;

  action_collection actions_hirep;
  std::vector<ActPtr>& actions;

  explicit MyActionLevel(unsigned int mul = 1) : 
    actions(std::get<0>(actions_hirep)), multiplier(mul) {
    // initialize the hirep vectors to zero.
    // apply(this->resize, actions_hirep, 0); //need a working resize
    assert(mul >= 1);
  }

  template < class GenField >
  void push_back(MyAction<GenField>* ptr) {
    // insert only in the correct vector
    std::get< Index < GenField, action_hirep_types>::value >(actions_hirep).push_back(ptr);
  }

  template <class ActPtr>
  static void resize(ActPtr ap, unsigned int n) {
    ap->resize(n);
  }

  // Loop on tuple for a callable function
  template <std::size_t I = 1, typename Callable, typename ...Args>
  inline typename std::enable_if<I == std::tuple_size<action_collection>::value, void>::type apply(Callable, Repr& R,Args&...) const {}

  template <std::size_t I = 1, typename Callable, typename ...Args>
  inline typename std::enable_if<I < std::tuple_size<action_collection>::value, void>::type apply(Callable fn, Repr& R, Args&... arguments) const {
    fn(std::get<I>(actions_hirep), std::get<I>(R.rep), arguments...);
    apply<I + 1>(fn, R, arguments...);
  }  

};

// Define the MyActionSet
template <class GaugeField, class R>
using MyActionSet = std::vector<MyActionLevel<GaugeField, R> >;




//////////////////////////////////////// Actions ///////////////////////////////////


template <class Gimpl>
class My_WilsonAction : public MyAction<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  explicit My_WilsonAction(RealD beta_) : beta(beta_){};

  virtual std::string action_name() {return "WilsonGaugeAction";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[WilsonGaugeAction] Beta: " << beta << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  virtual RealD S(const GaugeField &U) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    return Waction.S(U);
  }

  virtual void deriv(const GaugeField &U, GaugeField &dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    Waction.deriv(U, dSdU);
  }
private:
  RealD beta;
};


template <class Gimpl>
class My_WilsonAction_cell : public MyAction<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  explicit My_WilsonAction_cell(RealD beta_) : beta(beta_){};

  virtual std::string action_name() {return "WilsonGaugeAction_cell";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[WilsonGaugeAction_cell] Beta: " << beta << std::endl;
    return sstream.str();
  }

  virtual void refresh(const GaugeField &U,
                       GridParallelRNG &pRNG){};  // noop as no pseudoferms

  virtual RealD S(const GaugeField &U) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    return Waction.S(U);

  }

  virtual void deriv(const GaugeField &U, GaugeField &dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) {
    // std::cout << "before deriv" << std::endl;
    WilsonGaugeAction<Gimpl> Waction(beta);
    Waction.deriv(U, dSdU);

    // LatticeLorentzScalar mask = get_mask(U.Grid(), cell_size);
    if(mask == NULL) {
      mask = new LatticeLorentzScalar(U.Grid());
      Coordinate cell_size = pRNG.Grid()->_fdimensions;
      *mask = get_mask(U.Grid(), cell_size);
    }
    dSdU = dSdU * (*mask);
    // std::cout << "after deriv" << std::endl;
  }
private:
  RealD beta;
  LatticeLorentzScalar *mask=NULL;
};







template <class Gimpl>
class GFAction : public MyAction<typename Gimpl::GaugeField> {
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

  virtual void deriv(const GaugeField &U, GaugeField &dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) {
    WilsonGaugeAction<Gimpl> Waction(beta);
    GaugeField dSwdU(U.Grid());
    Waction.deriv(U, dSwdU);

  	RealD factor = 0.5 * (1./3.) * betaMM;

    GaugeField dSGF1dU(U.Grid());
    dSGF1dU = factor * Ta(U);

    // dSdU = dSwdU + dSGF1dU;
    // std::cout << "not adding dSGF2dU" << std::endl;

    GaugeField dSGF2dU(U.Grid());
    dSGF2dU = Zero();
    static LatticeColourMatrix g(U.Grid());
    if(first_step) {    // If it is the first step in a trajectory, set g=1.0 and run extra 10 heatbath sweeps
      g = 1.0;
      int sweeps = 10;
      GF_heatbath(U, g, sweeps, betaMM, table_path, pRNG); 
    }

    GF_heatbath(U, g, hb_offset, betaMM, table_path, pRNG); //hb_nsweeps before calculate equilibrium value
    GF_heatbath(U, g, innerMC_N, betaMM, table_path, pRNG, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU

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
class GFAction_cell : public MyAction<typename Gimpl::GaugeField> {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  explicit GFAction_cell(RealD beta_, RealD betaMM_, int innerMC_N_, int hb_offset_, const std::string& table_path_)
  : beta(beta_), betaMM(betaMM_), innerMC_N(innerMC_N_), hb_offset(hb_offset_), table_path(table_path_){}

  virtual std::string action_name() {return "GFAction_cell";}

  virtual std::string LogParameters(){
    std::stringstream sstream;
    sstream << GridLogMessage << "[GFAction_cell] Beta: " << beta  <<  "\t BetaMM: " << betaMM
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

  virtual void deriv(const GaugeField &U, GaugeField &dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG_cell, bool first_step) {
    GridBase *cell_grid = pRNG_cell.Grid();
    Coordinate cell_size = cell_grid->_fdimensions;

    WilsonGaugeAction<Gimpl> Waction(beta);
    GaugeField dSwdU(U.Grid());
    Waction.deriv(U, dSwdU);

    RealD factor = 0.5 * (1./3.) * betaMM;

    GaugeField dSGF1dU(U.Grid());
    dSGF1dU = factor * Ta(U);

    //////// Extract U_cell
    LatticeGaugeField U_cell(cell_grid);
    localCopyRegion(U, U_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
    // LatticeLorentzScalar cell_mask = get_cell_mask(cell_grid);
    // U_cell = U_cell * cell_mask;
    if(cell_mask == NULL) {
      cell_mask = new LatticeLorentzScalar(cell_grid);
      *cell_mask = get_cell_mask(cell_grid);
    }
    U_cell = U_cell * (*cell_mask);

    // LatticeColourMatrix g_cell(cell_grid); g_cell = 1.0;
    static LatticeColourMatrix g_cell(cell_grid);
    if(first_step) {    // If it is the first step in a trajectory, set g=1.0 and run extra 10 heatbath sweeps
      g_cell = 1.0;
      int sweeps = 10;
      GF_heatbath(U_cell, g_cell, sweeps, betaMM, table_path, pRNG_cell); //hb_nsweeps before calculate equilibrium value
    }

    GaugeField dSGF2dU_cell(cell_grid); dSGF2dU_cell = Zero();
    GF_heatbath(U_cell, g_cell, hb_offset, betaMM, table_path, pRNG_cell); //hb_nsweeps before calculate equilibrium value
    GF_heatbath(U_cell, g_cell, innerMC_N, betaMM, table_path, pRNG_cell, &dSGF2dU_cell, dOmegadU_g); // calculate dSGF2dU

    GaugeField dSGF2dU(U.Grid()); dSGF2dU = Zero();
    localCopyRegion(dSGF2dU_cell, dSGF2dU, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

    dSGF2dU = factor *  (1.0 / double(innerMC_N)) * dSGF2dU;

    dSdU = dSwdU + dSGF1dU - dSGF2dU;

    // LatticeLorentzScalar mask = get_mask(U.Grid(), cell_size);
    // dSdU = dSdU * mask;
    if(mask == NULL) {
      mask = new LatticeLorentzScalar(U.Grid());
      *mask = get_mask(U.Grid(), cell_size);
    }
    dSdU = dSdU * (*mask);
  }
private:
  RealD beta;
  RealD betaMM;
  int innerMC_N;
  int hb_offset;
  std::string table_path;
  LatticeLorentzScalar *cell_mask=NULL;
  LatticeLorentzScalar *mask=NULL;
};



// template <class Gimpl>
// class GF_DBW2Action : public MyAction<typename Gimpl::GaugeField> {
//  public:
//   INHERIT_GIMPL_TYPES(Gimpl);
//
//   explicit GF_DBW2Action(RealD beta_, RealD betaMM_, int innerMC_N_, int hb_offset_, const std::string& table_path_)
//   : beta(beta_), betaMM(betaMM_), innerMC_N(innerMC_N_), hb_offset(hb_offset_), table_path(table_path_){}
//
//   virtual std::string action_name() {return "GF_DBW2Action";}
//
//   virtual std::string LogParameters(){
//     std::stringstream sstream;
//     sstream << GridLogMessage << "[GF_DBW2Action] Beta: " << beta  <<  "\t BetaMM: " << betaMM
//       << "\t innerMC_N: " <<innerMC_N << "\t hb_offset: " << hb_offset <<  std::endl;
//     return sstream.str();
//   }
//
//   virtual void refresh(const GaugeField &U,
//                        GridParallelRNG &pRNG){};  // noop as no pseudoferms
//
//   //delta S_GF2 is calculated in GF_HMC.h: evolve();
//   virtual RealD S(const GaugeField &U) {
//     DBW2GaugeAction<Gimpl> DBW2_action(beta);
//
//     RealD Sw = DBW2_action.S(U);
//   	std::cout  << "DBW2 S: " <<  std::setprecision(15) << Sw << std::endl;
//
//     RealD SGF1 = - betaMM * Omega_no_g(U);
//   	std::cout  << "SGF1: " <<  std::setprecision(15) << SGF1 << std::endl;
//
//     return Sw + SGF1;
//   }
//
//   // virtual void deriv(const GaugeField &U, GaugeField &dSdU) {
//   virtual void deriv(const GaugeField &U, GaugeField &dSdU, GridSerialRNG &sRNG, GridParallelRNG &pRNG) {
//
//     DBW2GaugeAction<Gimpl> DBW2_action(beta);
//     GaugeField dSwdU(U.Grid());
//     DBW2_action.deriv(U, dSwdU);
//
//   	RealD factor = 0.5 * betaMM;
//
//     GaugeField dSGF1dU(U.Grid());
//     dSGF1dU = factor * Ta(U);
//
//     GaugeField dSGF2dU(U.Grid());
//     dSGF2dU = Zero();
//     LatticeColourMatrix g(U.Grid());
//     g = 1.0;
//     GF_heatbath(U, g, hb_offset, betaMM, table_path, pRNG); //hb_nsweeps before calculate equilibrium value
//
//     GF_heatbath(U, g, innerMC_N, betaMM, table_path, pRNG, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU
//
//     dSGF2dU = factor *  (1.0 / double(innerMC_N)) * dSGF2dU;
//
//     dSdU = dSwdU + dSGF1dU - dSGF2dU;
//
//   }
// private:
//   RealD beta;
//   RealD betaMM;
//   int innerMC_N;
//   int hb_offset;
//   std::string table_path;
// };

typedef My_WilsonAction<PeriodicGimplR>          My_WilsonActionR;
typedef My_WilsonAction_cell<PeriodicGimplR>          My_WilsonActionR_cell;
typedef GFAction<PeriodicGimplR>          GFActionR;
typedef GFAction_cell<PeriodicGimplR>          GFActionR_cell;
// typedef GF_DBW2Action<PeriodicGimplR>     GF_DBW2ActionR;

}
