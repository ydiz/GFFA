#include "../GFFA.h"

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

  // HMC_PARA hmc_para {};  // initialize each data member to default value
  // init(argc, argv, hmc_para);

  JSONReader reader("GFFA.json");

  GFFAParams hmc_para(reader);
  std::cout << hmc_para << std::endl; // FIXME:


  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({8,8,8,8}), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
  // GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(Coordinate({16,16,16,16}), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  GridParallelRNG _pRNG(grid);  // not used

  _pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

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
  // if(argc < 1) {
  //   std::cout << "You have to use input argv[1] as base directory for configuration" << std::endl;
  // }

  double interval = 0.1;
  // bool display_each_force = false;
  bool display_each_force = true;

  // std::string base_dir = argv[1];
  std::string base_dir = "/home/ahmedsheta/cuth_runs/results/GFFA_runs/bad_polyakov_lines/beta10/M3/traj0.6/all_steps24_MC40";
  // std::string base_dir = ".";
  // int traj_start = 1990, traj_end = 1995, traj_sep = 1; // for 24ID, kaon wall
  // int traj_start = 1000, traj_end = 1005, traj_sep = 1; // for 24ID, kaon wall
  int traj_start = 1000, traj_end = 1000, traj_sep = 1; // for 24ID, kaon wall
  // std::string base_dir = "/home/ydzhao/cuth/GFFA/results/GFFA_runs/beta=100/M3_eps0.1";
  // std::string base_dir = "/home/ydzhao/cuth/GFFA/beta100/FA_eps0.1";
  // int traj_start = 1000, traj_end = 1000, traj_sep = 1; // for 24ID, kaon wall
  // int traj_start = 5000, traj_end = 5000, traj_sep = 10; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    LatticeGaugeField U(grid);
    readField(U, base_dir + "/ckpoint_lat." + std::to_string(traj));  // ./ckpoint_lat.1000: 16nt16 lattice
    // U = 1.0;          // FIXME: I am setting to cold configuration

    // U = Zero();  
    // autoView( U_v, U, CpuWrite);
    // accelerator_for(ss, grid->oSites(),1,
    // {
    //   for(int mu=0; mu<4; ++mu) {
    //   U_v[ss](mu)()(0, 1) = 1.0;
    //   U_v[ss](mu)()(1, 0) = -1.0;
    //   U_v[ss](mu)()(2, 2) = 1.0;
    //   }
    // });
    //
    // PeriodicGimplR::HotConfiguration(pRNG, U);
    // print_grid_field_site(U, {0,0,0,0});
    // print_grid_field_site(U, {1,0,0,0});

    // LatticeGaugeField force(grid);
    // action->deriv(U, force); // force contains coefficient

    WilsonGaugeAction<PeriodicGimplR> Waction(hmc_para.beta);
    LatticeGaugeField dSwdU(U.Grid());
    Waction.deriv(U, dSwdU);

  	RealD factor = 0.5 * hmc_para.betaMM;

    LatticeGaugeField dSGF1dU(U.Grid());
    dSGF1dU = factor * Ta(U);

    LatticeGaugeField dSGF2dU(U.Grid());
    dSGF2dU = Zero();
    static LatticeColourMatrix g(U.Grid());
    static bool g_initialized = false;
    if(! g_initialized) {
      g = 1.0;
      g_initialized = true;
    }

    GF_heatbath(U, g, hmc_para.hb_offset, hmc_para.betaMM, hmc_para.table_path, _pRNG); //hb_nsweeps before calculate equilibrium value
    GF_heatbath(U, g, hmc_para.innerMC_N, hmc_para.betaMM, hmc_para.table_path, _pRNG, &dSGF2dU, dOmegadU_g); // calculate dSGF2dU

    dSGF2dU = factor *  (1.0 / double(hmc_para.innerMC_N)) * dSGF2dU;


    // std::cout << "Force dSGF2dU" << std::endl;
    // print_grid_field_site(dSGF2dU, {0,0,0,0});
    // std::cout << dSGF2dU << std::endl;
    // return 0;


    // std::cout << "===================dSGF2dU force: =====================" << std::endl;
    // get_force_stats(dSGF2dU, interval, hmc_para.epsilon);

    LatticeGaugeField dSdU(U.Grid());
    dSdU = dSwdU + dSGF1dU - dSGF2dU;

    if(display_each_force) {
      std::cout << "===================Wilson force: =====================" << std::endl;
      get_force_stats(dSwdU, interval, hmc_para.epsilon);

      std::cout << "===================dSGF1dU force: =====================" << std::endl;
      get_force_stats(dSGF1dU, interval, hmc_para.epsilon);

      std::cout << "===================dSGF2dU force: =====================" << std::endl;
      get_force_stats(dSGF2dU, interval, hmc_para.epsilon);
    }

    std::cout << "===================total force: =====================" << std::endl;
    get_force_stats(dSdU, interval, hmc_para.epsilon);
    //
    // // FFT theFFT((Grid::GridCartesian *)grid);
    // // theFFT.FFT_all_dim(force, force, FFT::forward);
    // //
    // //
    // // LatticeGaugeField force_L(grid);
    // // LatticeGaugeField force_T(grid);
    // // force_L = PL_projection(force, hmc_para.epsilon);
    // // force_T = force - force_L;
    // //
    // // std::cout << "force avg: " << get_average_force(force, interval) << std::endl;
    // // std::cout << "force_L avg: " << get_average_force(force_L, interval) << std::endl;
    // // std::cout << "force_T avg: " << get_average_force(force_T, interval) << std::endl;
  }

  // LatticeColourMatrix tmp(grid);
  // norm2(tmp);

  std::cout << "Finished!"  << std::endl;
  Grid_finalize();

} // main
