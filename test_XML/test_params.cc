#include <Grid/Grid.h>
// #include <iostream>
// #include <sys/sysinfo.h>


#include "../GF_params.h"
// #include "GF_Util.h"
// #include "GF_init_k.h"
// #include "measure_A.h"
#include "../observable.h"
// #include "GF_assert.h"
// #include "Integral_table.h"
// #include "subgroup_hb.h"
// #include "GF_heatbath_Util.h"
// #include "GF_generate_P.h"
// #include "GF_Action.h"
// #include "GF_deltaU.h"
// #include "GF_hmc_integrator.h"
// #include "GF_hmc_integrator_alg.h"
// #include "GF_HMC.h"
// #include "GF_GenericHMCrunner.h"


// #include <fenv.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char **argv) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);

  Grid_init(&argc, &argv);
  GridLogLayout();

  // JSONReader reader("GFFA.json");
  XmlReader reader("GFFA.xml");

  HMCparameters hmcParams(reader);
  GFFAParams hmc_para(reader);
  MyTC_para tc_para(reader);
  GaugeModes_para gm_para(reader);


  std::cout << gm_para << std::endl;

  // XmlWriter WRx("GFFA.xml");
  // write(WRx, "HMC", hmcParams);
  // write(WRx, "GFFA", hmc_para);
  // write(WRx, "GaugeModes_Observable", gm_para);
  // write(WRx, "WilsonFlow_Observable", tc_para);


  Grid_finalize();

} // main
