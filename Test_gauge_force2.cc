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

namespace Grid {
namespace QCD {
void localIndexToLocalGlobalCoor(GridBase *grid, int ss, std::vector<int> &lcoor, std::vector<int> &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  std::vector<int> processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}

LatticeGaugeField zyd_RandomGaugeTransform(GridParallelRNG &pRNG, const LatticeGaugeField &Umu){ // g = e^{i * scale \sum_a rv *  ta} // rv is U(-0.5, 0.5) rv
  double scale = 10.;
  LatticeColourMatrix g(Umu._grid);
  SU<3>::LieRandomize(pRNG, g, 5.0);

  LatticeGaugeField Utrans(Umu._grid);
  Utrans = Umu;
  SU<3>::GaugeTransform(Utrans, g);
  return Utrans;
}

}}

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
  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid({16,16,16,16}, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  GridParallelRNG   pRNG(grid); 
  pRNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));


  std::string base_dir = argv[1];
  int traj_start = 5000, traj_end = 5010, traj_sep = 10; // for 24ID, kaon wall
  // int traj_start = 5000, traj_end = 5000, traj_sep = 10; // for 24ID, kaon wall
  int traj_num = (traj_end - traj_start) / traj_sep + 1;


  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
    LatticeGaugeField U(grid);
    readField(U, base_dir + "/ckpoint_lat." + std::to_string(traj));


    LatticeGaugeField one(grid);
    one = 1.;
    
    for(int dummy = 0; dummy < 10; ++dummy) {
      LatticeGaugeField Utrans(grid);
      Utrans = zyd_RandomGaugeTransform(pRNG, U);
      Utrans = ProjectOnGroup(Utrans);  //Utrans is not exactly unitary, and gauge fixing will fail

      double alpha = 0.1;
      FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Utrans, alpha, 10000, 1.0e-8, 1.0e-10, true); // use FourierAccelSteepestDescentStep 

      LatticeGaugeField diff = Utrans - one;

      std::vector<double> vec(grid->lSites() * 4);
      parallel_for(int ss=0; ss<U._grid->lSites(); ss++) {
        std::vector<int> lcoor, gcoor;
        localIndexToLocalGlobalCoor(U._grid, ss, lcoor, gcoor);

        typename LatticeGaugeField::vector_object::scalar_object m;
        peekLocalSite(m, diff, lcoor);

        for(int mu=0; mu<4; ++mu) {
          typename LatticeColourMatrix::vector_object::scalar_object mm;
          mm() = m(mu);
          double val = norm2(mm);
          vec[ss * 4 + mu] = val;
        }
      }

      std::cout << vec << std::endl;
    }
  }
  // for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {
  //   LatticeGaugeField U(grid);
  //   readField(U, base_dir + "/ckpoint_lat." + std::to_string(traj));
  //
  //   // double alpha = 0.1;
  //   // // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U, alpha, 10000, 1.0e-12, 1.0e-12, false);
  //   // // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U, alpha, 10000, 1.0e-10, 1.0e-8, false);
  //   // FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U, alpha, 10000, 1.0e-10, 1.0e-8, true); // use FourierAccelSteepestDescentStep 
  //
  //   LatticeGaugeField one(grid);
  //   one = 1.;
  //   
  //   for(int dummy = 0; dummy < 10; ++dummy) {
  //     zyd_RandomGaugeTransform(pRNG, U);
  //
  //     LatticeGaugeField diff = U - one;
  //
  //     std::vector<double> vec(grid->lSites() * 4);
  //     parallel_for(int ss=0; ss<U._grid->lSites(); ss++) {
  //       std::vector<int> lcoor, gcoor;
  //       localIndexToLocalGlobalCoor(U._grid, ss, lcoor, gcoor);
  //
  //       typename LatticeGaugeField::vector_object::scalar_object m;
  //       peekLocalSite(m, diff, lcoor);
  //
  //       for(int mu=0; mu<4; ++mu) {
  //         typename LatticeColourMatrix::vector_object::scalar_object mm;
  //         mm() = m(mu);
  //         double val = norm2(mm);
  //         vec[ss * 4 + mu] = val;
  //       }
  //     }
  //
  //     std::cout << vec << std::endl;
  //   }
  // }

  Grid_finalize();

} // main
