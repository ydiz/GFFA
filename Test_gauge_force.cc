#include <Grid/Grid.h>
#include <sys/sysinfo.h>

#include "GF_Util.h"
#include "observable.h"

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



// #include <fenv.h>

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
  //
  //
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid({16,16,16,16}, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  LatticeGaugeField U(grid);
  readField(U, argv[1]);

  // action
  GFActionR GF_Wilson_action(hmc_para.beta, hmc_para.betaMM, hmc_para.innerMC_N, hmc_para.hb_offset, hmc_para.table_path);

  LatticeGaugeField force(grid);
  GF_Wilson_action.deriv(U, force); // force contains coefficient

  FFT theFFT((Grid::GridCartesian *)grid);
  theFFT.FFT_all_dim(force, force, FFT::forward);

  LatticeColourMatrix tmp(grid);
  norm2(tmp);

  Grid_finalize();

} // main
