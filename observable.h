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



class MyTC_para {
public:
  std::string type;
  double step_size;
  double adaptiveErrorTolerance;
  double maxTau;

  int TrajectoryStart;
  int TrajectoryInterval;

  bool saveSmearField;
  std::string smearFieldFilePrefix;
  std::string topoChargeOutFile;
};

std::ostream& operator<<(std::ostream &out, const MyTC_para &p) {
  out << "Topological Charge: "<< std::endl;
  out << "type: " << p.type << std::endl;
  out << "TrajectoryStart: " << p.TrajectoryStart << std::endl;
  out << "TrajectoryInterval: " << p.TrajectoryInterval << std::endl;
  out << "step_size: " << p.step_size << std::endl;
  out << "adaptiveErrorTolerance: " << p.adaptiveErrorTolerance << std::endl;
  out << "maxTau: " << p.maxTau << std::endl;
  out << "topoChargeOutFile: " << p.topoChargeOutFile << std::endl;
  out << "saveSmearField: " << p.saveSmearField << std::endl;
  out << "smearFieldFilePrefix: " << p.smearFieldFilePrefix << std::endl;
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

    MyWilsonFlow<PeriodicGimplR> WF(Par.step_size, Par.adaptiveErrorTolerance, Par.maxTau);

    if(traj > Par.TrajectoryStart && traj % Par.TrajectoryInterval == 0)
    {
      LatticeGaugeField Uflow(U._grid);

      if(Par.type=="fixedMaxTau") WF.smear_adaptive_fixed_tau(Uflow, U);
      else if(Par.type=="tSquaredE0.3") WF.smear_adaptive(Uflow, U);
      else assert(0);

      if(Par.saveSmearField) writeField(Uflow, Par.smearFieldFilePrefix + "." + std::to_string(traj));

      std::vector<double> topoCharge = timeSliceTopologicalCharge(Uflow);
      writeVector(topoCharge, traj, Par.topoChargeOutFile, U._grid->ThisRank());

      int def_prec = std::cout.precision();
      std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
        << "TC: [ " << traj << " ] : ";
      for(double x: topoCharge) std::cout << x << " ";
      std::cout << std::endl;

      std::cout.precision(def_prec);

    }
  }
};


// template < class Impl >
// class MyTCMod: public ObservableModule<MyTC<Impl>, MyTC_para>{
//   typedef ObservableModule<MyTC<Impl>, MyTC_para> ObsBase;
//   using ObsBase::ObsBase; // for constructors
//
//   // acquire resource
//   virtual void initialize(){
//     this->ObservablePtr.reset(new MyTC<Impl>());
//   }
//   public:
//   MyTCMod(): ObsBase(NoParameters()){}
// };


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
