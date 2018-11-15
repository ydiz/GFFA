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

  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  ActionLevel<HMCWrapper::Field> Level1(1);

  GFActionR Gaction(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.hb_multi_hit);
  WilsonGaugeActionR Waction(hmc_para.beta);

  //cannot define action inside if statement
  if(hmc_para.newAction){
    Level1.push_back(&Gaction);
  }
  else{
    Level1.push_back(&Waction);
  }

  TheHMC.TheAction.push_back(Level1);

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
