#include <Grid/Grid.h>
#include <sys/sysinfo.h>

#include "GF_Util.h"
// #include "observable.h"

#include "GF_HMC_para.h"
#include "GF_init.h"
#include "GF_assert.h"
#include "GF_init_k.h"
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

#include "Test_gauge_force.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  using namespace Grid;
  using namespace Grid::QCD;

  std::cout << "before Grid init" << std::endl;
  Grid_init(&argc, &argv);
  std::cout << "after Grid init" << std::endl;
  GridLogLayout();

  HMC_PARA hmc_para {};  // initialize each data member to default value
  init(argc, argv, hmc_para);


  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid({16,16,16,16}, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  // // action
  // Action<PeriodicGimplR::GaugeField> *action;
  // WilsonGaugeActionR Wilson_action(hmc_para.beta);
  // GFActionR GF_Wilson_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);
  // if(hmc_para.action == "Wilson"){
  //   action = &Wilson_action;
  // }
  // else if(hmc_para.action == "GF_Wilson"){
  //   action = &GF_Wilson_action;
  // }
  // else {
  //   std::cout << "Action not available" << std::endl;
  //   return 0;
  // }
  //
  if(argc < 1) {
    std::cout << "You have to use input argv[1] as base directory for configuration" << std::endl;
  }

  double interval = 0.1;

  std::string base_dir = argv[1];
  int traj_start = 5000, traj_end = 5010, traj_sep = 10; // for 24ID, kaon wall
  // int traj_start = 5000, traj_end = 5000, traj_sep = 10; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    LatticeGaugeField U(grid);
    readField(U, base_dir + "/ckpoint_lat." + std::to_string(traj));

    // LatticeGaugeField force(grid);
    // action->deriv(U, force); // force contains coefficient

    WilsonGaugeAction<PeriodicGimplR> Waction(hmc_para.beta);
    LatticeGaugeField dSwdU(U._grid);
    Waction.deriv(U, dSwdU);

  	RealD factor = 0.5 * hmc_para.betaMM;

    LatticeGaugeField dSGF1dU(U._grid);
    dSGF1dU = factor * Ta(U);

    LatticeGaugeField dSGF2dU(U._grid);
    dSGF2dU = zero;
    static LatticeColourMatrix g(U._grid);
    static bool g_initialized = false;
    if(! g_initialized) {
      g = 1.0;
      g_initialized = true;
    }

    GF_heatbath(U, g, hmc_para.hb_offset, hmc_para.betaMM, hmc_para.table_path); //hb_nsweeps before calculate equilibrium value
    GF_heatbath(U, g, hmc_para.innerMC_N, hmc_para.betaMM, hmc_para.table_path, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU

    dSGF2dU = factor *  (1.0 / double(hmc_para.innerMC_N)) * dSGF2dU;

    LatticeGaugeField dSdU(U._grid);
    dSdU = dSwdU + dSGF1dU - dSGF2dU;


    std::cout << "===================Wilson force: =====================" << std::endl;
    get_force_stats(dSwdU, interval, hmc_para.epsilon);

    std::cout << "===================dSGF1dU force: =====================" << std::endl;
    get_force_stats(dSGF1dU, interval, hmc_para.epsilon);

    std::cout << "===================dSGF2dU force: =====================" << std::endl;
    get_force_stats(dSGF2dU, interval, hmc_para.epsilon);

    std::cout << "===================total force: =====================" << std::endl;
    get_force_stats(dSdU, interval, hmc_para.epsilon);

    // FFT theFFT((Grid::GridCartesian *)grid);
    // theFFT.FFT_all_dim(force, force, FFT::forward);
    //
    //
    // LatticeGaugeField force_L(grid);
    // LatticeGaugeField force_T(grid);
    // force_L = PL_projection(force, hmc_para.epsilon);
    // force_T = force - force_L;
    //
    // std::cout << "force avg: " << get_average_force(force, interval) << std::endl;
    // std::cout << "force_L avg: " << get_average_force(force_L, interval) << std::endl;
    // std::cout << "force_T avg: " << get_average_force(force_T, interval) << std::endl;
  }

  // LatticeColourMatrix tmp(grid);
  // norm2(tmp);

  Grid_finalize();

} // main
