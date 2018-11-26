#include <Grid/Grid.h>
#include <iostream>
#include <sys/sysinfo.h>
#include "GF_HMC_para.h"
#include "GF_init.h"
// #include "GF"
#include "GF_assert.h"
#include "GF_Util.h"
#include "GF_init_k.h"
#include "GF_heatbath_Util.h"
#include "GF_generate_P.h"
#include "GF_Action.h"
#include "GF_deltaU.h"
#include "GF_hmc_integrator.h"
#include "GF_hmc_integrator_alg.h"
#include "GF_HMC.h"
#include "GF_GenericHMCrunner.h"

// #include <fenv.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

  HMC_PARA hmc_para {};  // initialize each data member to default value
  init(argc, argv, hmc_para);

  // typedef GF_GenericHMCRunner<GFMinimumNorm2> HMCWrapper;
  typedef GF_GenericHMCRunner<GFLeapFrog> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;

  TheHMC.Resources.AddFourDimGrid("gauge");

  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = hmc_para.saveInterval;
  CPparams.format = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TopologyObsParameters TopParams;
  TopParams.interval = hmc_para.TC_interval;
  TopParams.do_smearing = hmc_para.TC_do_smearing;
  TopParams.Smearing.steps = hmc_para.TC_Smearing_steps;
  TopParams.Smearing.step_size = hmc_para.TC_Smearing_step_size;
  TopParams.Smearing.meas_interval = hmc_para.TC_Smearing_meas_interval;
  TopParams.Smearing.maxTau = hmc_para.TC_Smearing_maxTau;
  // TopParams.interval = 5;
  // TopParams.do_smearing = true;
  // TopParams.Smearing.steps = 200;
  // TopParams.Smearing.step_size = 0.01;
  // TopParams.Smearing.meas_interval = 50;
  // TopParams.Smearing.maxTau = 2.0;
  TheHMC.Resources.AddObservable<QObs>(TopParams);


  // action
  ActionLevel<HMCWrapper::Field> Level1(1);

  WilsonGaugeActionR Wilson_action(hmc_para.beta);
  DBW2GaugeAction<PeriodicGimplR> DBW2_action(hmc_para.beta);
  GFActionR GF_Wilson_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.hb_multi_hit);

  if(hmc_para.action == "Wilson"){
    Level1.push_back(&Wilson_action);
  }
  else if(hmc_para.action == "GF_Wilson"){
    Level1.push_back(&GF_Wilson_action);
  }
  else if(hmc_para.action == "DBW2"){
    Level1.push_back(&DBW2_action);
  }
  else {
    std::cout << "Action not available" << std::endl;
    return 0;
  }

  TheHMC.TheAction.push_back(Level1);

  // HMC
  TheHMC.Parameters.NoMetropolisUntil = hmc_para.Thermalizations;
  TheHMC.Parameters.Trajectories = hmc_para.Trajectories;
  TheHMC.Parameters.MD.MDsteps = hmc_para.mdSteps;
  TheHMC.Parameters.MD.trajL   = hmc_para.trajL;
  TheHMC.Parameters.StartingType = hmc_para.StartingType;
  TheHMC.Parameters.StartTrajectory = hmc_para.StartingTrajectory;

  // TheHMC.ReadCommandLine(argc, argv);

  cout << hmc_para << endl;

  TheHMC.Run(hmc_para);

  Grid_finalize();

} // main
