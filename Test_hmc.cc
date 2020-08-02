#include <Grid/Grid.h>
#include <iostream>
#include <sys/sysinfo.h>

#include "GF_Util.h"
#include "GF_init_k.h"

#include "Test_A.h"

#include "observable.h"
#include "GF_HMC_para.h"
#include "GF_init.h"
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



namespace Grid {

struct GFFAParams : Serializable {
public:
  std::vector<int> fdims; // lattice size

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AParams,
    int, traj,
    std::string, lat,
    std::string, config,
    int, nhits,
    std::string, ensemble_name,
    std::string, prefix,
    bool, doTimeDilution
  );

  template <class ReaderClass>
  A2AParams(Reader<ReaderClass>& reader) {
    read(reader, "A2A", *this);
    GridCmdOptionIntVector(lat, fdims);
    std::cout << "fdims: " << fdims << std::endl;
    std::cout << *this << std::endl;
  }

};

}



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
  typedef GF_GenericHMCRunner<GFLeapFrog> HMCWrapper;
  // typedef GF_GenericHMCRunner<GFForceGradient> HMCWrapper;
  HMCWrapper TheHMC;

  TheHMC.Resources.AddFourDimGrid("gauge");



  std::cout << "before reading json" << << std::endl;
  JSONReader reader("GFFA.json");
  std::cout << "After reading json" << std::endl;

  TheHMC.Parameters.initialize(reader);
  std::cout << TheHMC.Parameters << std::endl;

  // Checkpointer definition
  // CheckpointerParameters CPparams;
  CheckpointerParameters CPparams(reader);
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  // CPparams.saveInterval = hmc_para.saveInterval;
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
  typedef MyTCMod<HMCWrapper::ImplPolicy> QObs; 
  TheHMC.Resources.AddObservable<QObs>(hmc_para.tc_para);
  if(hmc_para.isGFFA) {
    typedef GaugeModesMod<HMCWrapper::ImplPolicy> GMObs; 
    TheHMC.Resources.AddObservable<GMObs>(hmc_para.gm_para);
  }


  // action
  ActionLevel<HMCWrapper::Field> Level1(1);

  WilsonGaugeActionR Wilson_action(hmc_para.beta);
  GFActionR GF_Wilson_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);
  DBW2GaugeAction<PeriodicGimplR> DBW2_action(hmc_para.beta);
  GF_DBW2ActionR GF_DBW2_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);

  if(hmc_para.action == "Wilson"){
    Level1.push_back(&Wilson_action);
  }
  else if(hmc_para.action == "GF_Wilson"){
    Level1.push_back(&GF_Wilson_action);
  }
  else if(hmc_para.action == "DBW2"){
    Level1.push_back(&DBW2_action);
  }
  else if(hmc_para.action == "GF_DBW2"){
    Level1.push_back(&GF_DBW2_action);
  }
  else {
    std::cout << "Action not available" << std::endl;
    return 0;
  }

  TheHMC.TheAction.push_back(Level1);


  // // HMC
  // TheHMC.Parameters.NoMetropolisUntil = hmc_para.Thermalizations;
  // TheHMC.Parameters.Trajectories = hmc_para.Trajectories;
  // TheHMC.Parameters.MD.MDsteps = hmc_para.mdSteps;
  // TheHMC.Parameters.MD.trajL   = hmc_para.trajL;
  // TheHMC.Parameters.StartingType = hmc_para.StartingType;
  // TheHMC.Parameters.StartTrajectory = hmc_para.StartingTrajectory;

  // TheHMC.ReadCommandLine(argc, argv);

  cout << hmc_para << endl;

  TheHMC.Run(hmc_para);

  Grid_finalize();

} // main
