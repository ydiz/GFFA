#include <Grid/Grid.h>
#include <iostream>
#include <sys/sysinfo.h>


#include "GF_params.h"
#include "GF_Util.h"
#include "GF_init_k.h"
#include "measure_A.h"
#include "observable.h"
#include "GF_assert.h"
#include "Integral_table.h"
#include "subgroup_hb.h"
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

  Grid_init(&argc, &argv);
  GridLogLayout();

  JSONReader reader("GFFA.json");

  GFFAParams hmc_para(reader);
  HMCparameters hmcParams(reader);
  MyTC_para tc_para(reader);
  GaugeModes_para gm_para(reader);
  PolyakovLoop_para pl_para(reader);

  std::cout << hmcParams << std::endl;
  std::cout << hmc_para << std::endl;
  if(hmc_para.measure_A) {
    hmcParams.StartTrajectory = 0;
    hmcParams.NoMetropolisUntil = 1;
    hmcParams.Trajectories = 0;
    tc_para.TrajectoryStart = 999999;
  }

  typedef GF_GenericHMCRunner<GFLeapFrog> HMCWrapper;
  // typedef GF_GenericHMCRunner<GFMinimumNorm2> HMCWrapper;
  // typedef GF_GenericHMCRunner<GFForceGradient> HMCWrapper;

  HMCWrapper TheHMC(hmcParams);

  TheHMC.Resources.AddFourDimGrid("gauge"); // This is done in initialize

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

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  typedef LinkTraceMod<HMCWrapper::ImplPolicy> LTObs;
  TheHMC.Resources.AddObservable<LTObs>();
  typedef PolyakovLoopMod<HMCWrapper::ImplPolicy> PolyakovObs;
  TheHMC.Resources.AddObservable<PolyakovObs>(pl_para);

  typedef MyTCMod<HMCWrapper::ImplPolicy> QObs; 
  TheHMC.Resources.AddObservable<QObs>(tc_para);
  if(hmc_para.isGFFA) {
    typedef GaugeModesMod<HMCWrapper::ImplPolicy> GMObs; 
    TheHMC.Resources.AddObservable<GMObs>(gm_para);
  }


  // action
  MyActionLevel<HMCWrapper::Field> Level1(1);

  My_WilsonActionR Wilson_action(hmc_para.beta);
  // WilsonGaugeActionR Wilson_action(hmc_para.beta);
  GFActionR GF_Wilson_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);
  // DBW2GaugeAction<PeriodicGimplR> DBW2_action(hmc_para.beta);
  GF_DBW2ActionR GF_DBW2_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);

  if(hmc_para.action == "Wilson"){
    Level1.push_back(&Wilson_action);
  }
  else if(hmc_para.action == "GF_Wilson"){
    Level1.push_back(&GF_Wilson_action);
  }
  // else if(hmc_para.action == "DBW2"){
  //   Level1.push_back(&DBW2_action);
  // }
  else if(hmc_para.action == "GF_DBW2"){
    Level1.push_back(&GF_DBW2_action);
  }
  else {
    std::cout << "Action not available" << std::endl;
    return 0;
  }

  TheHMC.TheAction.push_back(Level1);

  TheHMC.Run(hmc_para);

  Grid_finalize();

} // main
