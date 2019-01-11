#include <cmath>

namespace Grid {
namespace QCD {

template <class IntegratorType>
class GF_HybridMonteCarlo {
 private:
  const HMCparameters Params;

  typedef typename IntegratorType::Field Field;
  typedef std::vector< HmcObservable<Field> * > ObsListType;

  	//pass these from the resource manager
  GridSerialRNG &sRNG;
  GridParallelRNG &pRNG;

  Field &Ucur;

  IntegratorType &TheIntegrator;
	ObsListType Observables;

  bool metropolis_test(const RealD DeltaH) {
    RealD rn_test;

    RealD prob = std::exp(-DeltaH);

    random(sRNG, rn_test);

    std::cout << GridLogMessage
              << "--------------------------------------------------\n";
    std::cout << GridLogMessage << "exp(-dH) = " << prob
              << "  Random = " << rn_test << "\n";
    std::cout << GridLogMessage
              << "Acc. Probability = " << ((prob < 1.0) ? prob : 1.0) << "\n";

    if ((prob > 1.0) || (rn_test <= prob)) {  // accepted
      std::cout << GridLogMessage << "Metropolis_test -- ACCEPTED\n";
      std::cout << GridLogMessage
                << "--------------------------------------------------\n";
      return true;
    } else {  // rejected
      std::cout << GridLogMessage << "Metropolis_test -- REJECTED\n";
      std::cout << GridLogMessage
                << "--------------------------------------------------\n";
      return false;
    }
  }

  /////////////////////////////////////////////////////////
  // Evolution
  /////////////////////////////////////////////////////////
  RealD evolve_hmc_step(Field &U, const Momenta_k &KK) {
    TheIntegrator.GF_refresh(U, pRNG, KK);  // set U and initialize P and phi's

    // RealD H0 = TheIntegrator.GF_S(U, KK);  // initial state action
    //
    // std::streamsize current_precision = std::cout.precision();
    // std::cout.precision(15);
    // std::cout << GridLogMessage << "Total H before trajectory = " << H0 << "\n";
    // std::cout.precision(current_precision);

    TheIntegrator.GF_integrate(U, KK);

    // RealD H1 = TheIntegrator.GF_S(U, KK);  // updated state action
    //
    // std::cout.precision(15);
    // std::cout << GridLogMessage << "Total H after trajectory  = " << H1
	  //     << "  dH = " << H1 - H0 << "\n";
    // std::cout.precision(current_precision);

    // return (H1 - H0);
    return 0;
  }


 public:
  /////////////////////////////////////////
  // Constructor
  /////////////////////////////////////////
  GF_HybridMonteCarlo(HMCparameters _Pams, IntegratorType &_Int,
                   GridSerialRNG &_sRNG, GridParallelRNG &_pRNG,
                   ObsListType _Obs, Field &_U)
    : Params(_Pams), TheIntegrator(_Int), sRNG(_sRNG), pRNG(_pRNG), Observables(_Obs), Ucur(_U) {}
  ~GF_HybridMonteCarlo(){};

  void evolve(const HMC_PARA &HMC_para) {

    const Momenta_k KK(Ucur._grid, HMC_para.M, HMC_para.epsilon, HMC_para.newHp);

    Real DeltaH;

    Field Ucopy(Ucur._grid);

    Params.print_parameters();
    TheIntegrator.print_actions();

  	LatticeColourMatrix g(Ucur._grid);
  	g = 1.0;
    // Actual updates (evolve a copy Ucopy then copy back eventually)
    unsigned int FinalTrajectory = Params.Trajectories + Params.NoMetropolisUntil + Params.StartTrajectory;
    for (int traj = Params.StartTrajectory; traj < FinalTrajectory; ++traj) {
      std::cout << GridLogMessage << "-- # Trajectory = " << traj << "\n";
      if (traj < Params.StartTrajectory + Params.NoMetropolisUntil) {
      	std::cout << GridLogMessage << "-- Thermalization" << std::endl;
      }

      double t0=usecond();
      Ucopy = Ucur;

      DeltaH = evolve_hmc_step(Ucopy, KK);

      // no Metropolis, always accept.
      Ucur = Ucopy;

  	  // Real tt=0;
  	  // Real FinalDeltaH=0;

      //calculate delta S_GF2
      // if(HMC_para.newAction){
      //   //LatticeColourMatrix g(Ucur._grid);
      //   //g = 1.0;
      //   // GF_metro(Ucur, g, HMC_para.hb_offset, HMC_para.betaMM, HMC_para.hb_multi_hit, 0.05);
      //   GF_heatbath(Ucur, g, HMC_para.hb_offset, HMC_para.betaMM, HMC_para.hb_multi_hit);
      //   	std::cout << "Omega_g(g, Ucur): "<< Omega_g(g, Ucur) << std::endl;
      //   	std::cout << "Omega_g(g, Ucopy): "<< Omega_g(g, Ucopy) << std::endl;
      //   Real DeltaH_SG2=0;
      //   for(int i=0; i<HMC_para.innerMC_N; ++i)
      //   {
      //   	  tt = DeltaH + HMC_para.betaMM * (Omega_g(g, Ucopy) - Omega_g(g, Ucur));
      //   	  DeltaH_SG2 += std::exp( tt );
      //   	  if(i%100==0) std::cout << "DeltaH tt: "<< tt << std::endl;
      //     // GF_metro(Ucur, g, HMC_para.hb_nsweeps, HMC_para.betaMM, HMC_para.hb_multi_hit, 0.05);
      //     GF_heatbath(Ucur, g, HMC_para.hb_nsweeps, HMC_para.betaMM, HMC_para.hb_multi_hit);
      //   }
      //   FinalDeltaH = std::log(DeltaH_SG2/double(HMC_para.innerMC_N));
      //   //FinalDeltaH = DeltaH_SG2/double(HMC_para.innerMC_N);
      //   	std::cout<< GridLogMessage << "Final deltaH: " << FinalDeltaH << std::endl;
      // }

      // Metropolis-Hastings test
      // bool accept = true;
      // if (traj >= Params.StartTrajectory + Params.NoMetropolisUntil) {
      //   //accept = metropolis_test(DeltaH);
      //   accept = metropolis_test(FinalDeltaH);
      // } else {
      // 	std::cout << GridLogMessage << "Skipping Metropolis test" << std::endl;
      // }
      //
      // if (accept)
      //   Ucur = Ucopy;

      double t1=usecond();
      std::cout << GridLogMessage << "Total time for trajectory (s): " << (t1-t0)/1e6 << std::endl;

      for (int obs = 0; obs < Observables.size(); obs++) {
      	std::cout << GridLogDebug << "Observables # " << obs << std::endl;
      	std::cout << GridLogDebug << "Observables total " << Observables.size() << std::endl;
      	std::cout << GridLogDebug << "Observables pointer " << Observables[obs] << std::endl;
        Observables[obs]->TrajectoryComplete(traj + 1, Ucur, sRNG, pRNG);
      }
      std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    }
  }

};


}  // QCD
}  // Grid
