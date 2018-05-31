#include <Grid/Grid.h>
#include <iostream>
#include "GF_HMC_para.h"
#include "GF_init.h"
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

  int noMetro = 10;
  int traj = 21; //number of trajectories
  int mdSteps = 10;
  Real trajL = 1.0;

  HMC_PARA HMC_para;

  GF_init(argc, argv, noMetro, traj, mdSteps, trajL, HMC_para);

  typedef GF_GenericHMCRunner<GFMinimumNorm2> HMCWrapper;
  // typedef GF_GenericHMCRunner<GFLeapFrog> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;

  TheHMC.Resources.AddFourDimGrid("gauge");

  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 30;
  CPparams.format = "IEEE64BIG";

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  ActionLevel<HMCWrapper::Field> Level1(1);

  GFActionR Gaction(HMC_para.beta, HMC_para.betaMM, HMC_para.innerMC_N, HMC_para.hb_offset, HMC_para.hb_nsweeps);
  WilsonGaugeActionR Waction(HMC_para.beta);

  //cannot define action inside if statement
  if(HMC_para.newAction){
    Level1.push_back(&Gaction);
  }
  else{
    Level1.push_back(&Waction);
  }

  TheHMC.TheAction.push_back(Level1);

  TheHMC.Parameters.NoMetropolisUntil = noMetro;
  TheHMC.Parameters.Trajectories = traj;
  TheHMC.Parameters.MD.MDsteps = mdSteps;
  TheHMC.Parameters.MD.trajL   = trajL;
  // TheHMC.Parameters.StartingType = "ColdStart";

  TheHMC.ReadCommandLine(argc, argv);

  cout << HMC_para << endl;

  TheHMC.Run(HMC_para);

  Grid_finalize();

} // main
