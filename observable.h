#pragma once

#include "./WilsonFlow/WF.h"
#include "polyakov_util.h"

namespace Grid {

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




struct PolyakovLoop_para : Serializable {

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(PolyakovLoop_para,
    bool, do_measure
  );

  template <class ReaderClass>
  PolyakovLoop_para(Reader<ReaderClass>& reader) {
    read(reader, "PolyakovLoop_Observable", *this);
    // std::cout << *this << std::endl;
  }

};








template <class Impl>
class PolyakovLoopLogger : public HmcObservable<typename Impl::Field> {
  PolyakovLoop_para Par;
public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef typename Impl::Field Field;

  PolyakovLoopLogger(PolyakovLoop_para P): Par(P) {}

  void TrajectoryComplete(int traj,
                          Field &U,
                          GridSerialRNG &sRNG,
                          GridParallelRNG &pRNG) {

    if(!Par.do_measure) return;

    Coordinate fdims = U.Grid()->FullDimensions();
    int V = 1;
    for(int i=0; i<4; ++i) V *= fdims[i];
    // int T = fdims[3];

    Lattice<iVector<iScalar<iScalar<vComplex>>, 4>> all_polya(U.Grid());
    for(int mu=0; mu<4; ++mu) {       // Direction of Polaykov loop
      // LatticeColourMatrix tmp = peekLorentz(U, mu); // U_mu
      // for(int i=0; i<fdims[mu]-1; ++i) tmp = tmp * Cshift(tmp, mu, 1);  // WRONG
      LatticeColourMatrix Umu = peekLorentz(U, mu); // U_mu
      LatticeColourMatrix Umu_shift = Umu;
      LatticeColourMatrix tmp = Umu;
     
      for(int i=0; i<fdims[mu]-1; ++i) {
        Umu_shift = Cshift(Umu_shift, mu, 1); 
        tmp = tmp * Umu_shift;
      }   

      LatticeComplex polya_lat = trace(tmp);
      pokeLorentz(all_polya, polya_lat, mu);

        
      Complex avg_polya = sum(polya_lat)()()() / double(V);
      std::cout << "Average polyakov line for mu=[" << mu << "]: " << avg_polya << std::endl;

      // LatticeComplex polya_phase(U.Grid()), polya_norm(U.Grid());  // imaginary part is zero; 
      // lat_cmpl_to_arg_norm(polya_lat, polya_phase, polya_norm);
      //
      // double avg, stdev;
      //
      // avg_std_phases_norms(polya_phase, avg, stdev);
      // std::cout << "Phase avg for mu=[" << mu << "]: " << avg << std::endl;
      // std::cout << "Phase std for mu=[" << mu << "]: " << stdev << std::endl;
      // std::cout << "Phase bins for mu=[" << mu << "]: " << binning_phases(polya_phase) << std::endl;
      //
      // avg_std_phases_norms(polya_norm, avg, stdev);
      // std::cout << "Norm avg for mu=[" << mu << "]: " << avg << std::endl;
      // std::cout << "Norm std for mu=[" << mu << "]: " << stdev << std::endl;
      // std::cout << "Norm bins for mu=[" << mu << "]: " << binning_norms(polya_norm) << std::endl;

      // std::vector<Complex> polya(V);
      // autoView(polya_lat_v, polya_lat, CpuRead);
      // thread_for(ss, U.Grid()->lSites(), {
      //   Coordinate lcoor, gcoor;
      //   localIndexToLocalGlobalCoor(U.Grid(), ss, lcoor, gcoor);
      //
      //   if(gcoor[mu] == 0) {
      //     typename LatticeColourMatrix::vector_object::scalar_object m;
      //     peekLocalSite(m, tmp_v, lcoor);
      //
      //     typename LatticeComplex::vector_object::scalar_object m2;
      //     m2 = trace(m);
      //
      //     // int x = gcoor[0], y = gcoor[1], z = gcoor[2];
      //     // polya[x * fdims[0] * fdims[1] + y * fdims[1] + z] = m2()()();
      //
      //     int idx = gcoor[dir1] * fdims[dir1] * fdims[dir2] + gcoor[dir2] * fdims[dir2] + gcoor[dir3];
      //     polya[idx] = m2()()();
      //   }
      // });
      // int def_prec = std::cout.precision();
      // std::cout << std::setprecision(3) << "Polyakov loop: [ " << traj << " ] mu = " << mu << ": " << polya << std::endl;
      //
      // Complex avg_polya = 0.;
      // for(int i=0; i<polya.size(); ++i) avg_polya += polya[i];
      // avg_polya /= double(polya.size());
      //
      // std::cout << GridLogMessage
      //     << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
      //     << "AvgPolyakovLoop: [ " << traj << " ] mu = " << mu << ": " << avg_polya << std::endl;
      // std::cout.precision(def_prec);
    }
    writeScidac(all_polya, "polya."+std::to_string(traj));

  }
};

template < class Impl >
class PolyakovLoopMod: public ObservableModule<PolyakovLoopLogger<Impl>, PolyakovLoop_para>{
  typedef ObservableModule<PolyakovLoopLogger<Impl>, PolyakovLoop_para> ObsBase;
  using ObsBase::ObsBase; // for constructors

  // acquire resource
  virtual void initialize(){
    this->ObservablePtr.reset(new PolyakovLoopLogger<Impl>(this->Par_));
  }
public:
  PolyakovLoopMod(PolyakovLoop_para Par): ObsBase(Par){}
  // PolyakovLoopMod(): ObsBase(){}
};

// template < class Impl >
// class GaugeModesMod: public ObservableModule<GaugeModesLogger<Impl>, GaugeModes_para>{
//   typedef ObservableModule<GaugeModesLogger<Impl>, GaugeModes_para> ObsBase;
//   using ObsBase::ObsBase; // for constructors
//
//   // acquire resource
//   virtual void initialize() {
//     this->ObservablePtr.reset(new GaugeModesLogger<Impl>(this->Par_));
//   }
//   public:
//   GaugeModesMod(GaugeModes_para Par): ObsBase(Par){}
//   GaugeModesMod(): ObsBase(){}
// };













struct GaugeModes_para : Serializable {

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeModes_para,
    std::vector<std::vector<int>>, modes
  );

  template <class ReaderClass>
  GaugeModes_para(Reader<ReaderClass>& reader) {
    read(reader, "GaugeModes_Observable", *this);
    // std::cout << *this << std::endl;
  }

};



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

    using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
    static LatticeUeval last_log_evals(U.Grid()); 
    static bool last_log_evals_initialized = false;

    // measure_A(U, Par.coors, last_log_evals, last_log_evals_initialized);
    measure_A(U, Par.modes, last_log_evals, last_log_evals_initialized);

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
  // GaugeModesMod(): ObsBase(){}
};






struct MyTC_para: Serializable {

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MyTC_para,
    bool, do_measure,
    std::string, type,
    std::vector<double>, meas_taus,
    int, TrajectoryStart,
    int, TrajectoryInterval,
    double, step_size,
    double, adaptiveErrorTolerance
  );

  template <class ReaderClass>
  MyTC_para(Reader<ReaderClass>& reader) {
    read(reader, "WilsonFlow_Observable", *this);
    // std::cout << *this << std::endl;
  }

};





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

    if(!Par.do_measure) return;

    // MyWilsonFlow<PeriodicGimplR> WF(Par.step_size, Par.adaptiveErrorTolerance, Par.meas_taus, Par.topoChargeOutFile, traj);
    MyWilsonFlow<PeriodicGimplR> WF(Par.step_size, Par.adaptiveErrorTolerance, Par.meas_taus, traj);

    if(traj > Par.TrajectoryStart && traj % Par.TrajectoryInterval == 0)
    {
      LatticeGaugeField Uflow(U.Grid());

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
  // MyTCMod(): ObsBase(){}
};






}
