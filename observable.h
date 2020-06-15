#pragma once

#include "./WilsonFlow/WF.h"

namespace Grid {
namespace QCD {

template <class Impl>
class LinkTraceLogger : public HmcObservable<typename Impl::Field> {
public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef typename Impl::Field Field;

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    RealD lt = WilsonLoops<Impl>::linkTrace(U);

    int def_prec = std::cout.precision();

    std::cout << GridLogMessage
        << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "LinkTrace: [ " << traj << " ] "<< lt << std::endl;

    std::cout.precision(def_prec);

  }
};


template < class Impl >
class LinkTraceMod: public ObservableModule<LinkTraceLogger<Impl>, NoParameters>{
  typedef ObservableModule<LinkTraceLogger<Impl>, NoParameters> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new LinkTraceLogger<Impl>());
  }
public:
  LinkTraceMod(): ObsBase(NoParameters()){}
};











class GaugeModes_para {
public:
  std::vector<std::vector<int>> coors;
};

std::ostream& operator<<(std::ostream &out, const GaugeModes_para &p) {
  out << "GaugeModes: "<< std::endl;
  out << "coors: " << p.coors << std::endl;
  return out;
}




template <class Impl>
class GaugeModesLogger : public HmcObservable<typename Impl::Field> {
  GaugeModes_para Par;
public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef typename Impl::Field Field;

  GaugeModesLogger(GaugeModes_para P): Par(P) {}

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    // using LatticeGaugeFieldSite = typename LatticeGaugeField::vector_object::scalar_object;
    // using LatticeUevalSite = iVector<iScalar<iVector<vComplex, 3> >, 4>;
    // static std::vector<LatticeUevalSite> log_eval(coors.size(), 0.);
    //
    // static bool last_log_evals_initialized = false;

    using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
    static LatticeUeval last_log_evals(U._grid); 
    static bool last_log_evals_initialized = false;

    measure_A(U, Par.coors, last_log_evals, last_log_evals_initialized);

    // for(int i=0; i<Par.coors.size(); ++i) {
    //
    //   std::vector<int> coor = Par.coors[i];
    //   LatticeGaugeFieldSite m;
    //   peekSite(m, U, coor);
    //   
    //   for(int mu=0; mu<4; ++mu) {
    //     m(mu)() = Log( m(mu)(), log_eval[i](mu)(), last_log_evals_initialized);
    //   }
    //   last_log_evals_initialized = true;
    //   std::cout << GridLogMessage << "GaugeModes: [ " << traj << " ] coor [" << log_m << "]" << std::endl;
    // }
    //
    

  }
};


template < class Impl >
class GaugeModesMod: public ObservableModule<GaugeModesLogger<Impl>, GaugeModes_para>{
  typedef ObservableModule<GaugeModesLogger<Impl>, GaugeModes_para> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize() {
    this->ObservablePtr.reset(new GaugeModesLogger<Impl>(this->Par_));
  }
  public:
  GaugeModesMod(GaugeModes_para Par): ObsBase(Par){}
  GaugeModesMod(): ObsBase(){}
};






//
// template <class Impl>
// class WilsonLineTraceLogger : public HmcObservable<typename Impl::Field> {
// public:
//   INHERIT_GIMPL_TYPES(Impl);
//   typedef typename Impl::Field Field;
//
//   void TrajectoryComplete(int traj,
//                           Field &U,
//                           GridSerialRNG &sRNG,
//                           GridParallelRNG &pRNG) {
//
//
//     // !!! Cannot be parallel for
//     using LatticeGaugeFieldSite = typename LatticeGaugeField::vector_object::scalar_object;
//     int T = U._grid->_fdimensions[3];
//     std::vector<LatticeGaugeFieldSite> local_rst(T, 1.);
//     for(int ss=0; ss<force._grid->lSites(); ss++) {
//       std::vector<int> lcoor, gcoor;
//       localIndexToLocalGlobalCoor(U._grid, ss, lcoor, gcoor);
//       LatticeGaugeFieldSite m;
//       peekLocalSite(m, force, lcoor);
//       local_rst[gcoors[3]] *= m;
//     }
//
//     double rst;
//
//     int def_prec = std::cout.precision();
//
//     std::cout << GridLogMessage
//         << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
//         << "WilsonLineTrace: [ " << traj << " ] "<< rst<< std::endl;
//
//     std::cout.precision(def_prec);
//
//   }
// };
//
//
// template < class Impl >
// class WilsonLineTraceMod: public ObservableModule<WilsonLineTraceLogger<Impl>, NoParameters>{
//   typedef ObservableModule<WilsonLineTraceLogger<Impl>, NoParameters> ObsBase;
//   using ObsBase::ObsBase; // for constructors
//
//   // acquire resource
//   virtual void initialize(){
//     this->ObservablePtr.reset(new WilsonLineTraceLogger<Impl>());
//   }
// public:
//   WilsonLineTraceMod(): ObsBase(NoParameters()){}
// };









class MyTC_para {
public:
  std::string type;
  double step_size;
  double adaptiveErrorTolerance;
  // double maxTau;
  std::vector<double> meas_taus;

  int TrajectoryStart;
  int TrajectoryInterval;

  // bool saveSmearField;
  // std::string smearFieldFilePrefix;
  std::string topoChargeOutFile;
};

std::ostream& operator<<(std::ostream &out, const MyTC_para &p) {
  out << "Topological Charge: "<< std::endl;
  out << "type: " << p.type << std::endl;
  out << "TrajectoryStart: " << p.TrajectoryStart << std::endl;
  out << "TrajectoryInterval: " << p.TrajectoryInterval << std::endl;
  out << "step_size: " << p.step_size << std::endl;
  out << "adaptiveErrorTolerance: " << p.adaptiveErrorTolerance << std::endl;
  out << "meas_taus: " << p.meas_taus << std::endl;
  out << "topoChargeOutFile: " << p.topoChargeOutFile << std::endl;
  return out;
}


template <class Impl>
class MyTC : public HmcObservable<typename Impl::Field> {
  MyTC_para Par;
public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef typename Impl::Field Field;

  MyTC(MyTC_para P): Par(P) {}

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    MyWilsonFlow<PeriodicGimplR> WF(Par.step_size, Par.adaptiveErrorTolerance, Par.meas_taus, Par.topoChargeOutFile, traj);

    if(traj > Par.TrajectoryStart && traj % Par.TrajectoryInterval == 0)
    {
      LatticeGaugeField Uflow(U._grid);

      if(Par.type=="fixed_taus") WF.smear_adaptive_fixed_tau(Uflow, U);
      else if(Par.type=="tSquaredE0.3") WF.smear_adaptive_fixed0p3(Uflow, U);
      else assert(0);

    }
  }
};


template < class Impl >
class MyTCMod: public ObservableModule<MyTC<Impl>, MyTC_para>{
  typedef ObservableModule<MyTC<Impl>, MyTC_para> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new MyTC<Impl>(this->Par_));
  }
  public:
  MyTCMod(MyTC_para Par): ObsBase(Par){}
  MyTCMod(): ObsBase(){}
};






}}
