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

#include "Test_A.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

int main(int argc, char **argv) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  using namespace Grid;
  using namespace Grid::QCD;

  double beta = 10.;

  Grid_init(&argc, &argv);
  GridLogLayout();
  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid({16,16,16,16}, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  GridParallelRNG   pRNG(grid); 
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  // std::string base_dir = argv[1];
  std::string base_dir = "/home/ydzhao/cuth/gauge_force_config/beta10_M0.5_traj0.5";
  // int traj_start = 5000, traj_end = 5010, traj_sep = 10; // for 24ID, kaon wall
  int traj_start = 5000, traj_end = 5000, traj_sep = 10; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  // double interval = 1., epsilon = 0.2;
  double epsilon = 0.2;

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    LatticeGaugeField U(grid);
    readField(U, base_dir + "/ckpoint_lat." + std::to_string(traj));
    // measure_A(U, beta, interval, epsilon);
    measure_A(U, beta, epsilon);
    
    // LatticeGaugeField An(U._grid);
    // An = timesMinusI(Log(U)) * std::sqrt(6./beta);
    //
    // Lattice<iVector<iScalar<iVector<vComplex, 8 > >, 4> > alg(U._grid);
    //
    // for (int mu = 0; mu < Nd; mu++)
    // {
    //   LatticeColourMatrix An_mu(An._grid);
    //   SU3::LatticeAlgebraVector alg_mu(U._grid);
    //
    //   An_mu = peekLorentz(An, mu);
    //   SU3::projectOnAlgebra(alg_mu, An_mu);
    //   pokeLorentz(alg, alg_mu, mu);
    // }
    // // std::cout << alg << std::endl;
    //
    // for(int i=0; i<=8; ++i) {
    //   cout << "i=" << i << endl;
    //   print_grid_field_site(alg, {i,0,0,0});
    // }
    //

  }

  Grid_finalize();

} // main
