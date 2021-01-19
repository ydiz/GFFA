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

  GridParallelRNG pRNG(grid); 

  pRNG.SeedFixedIntegers(std::vector<int>({1,2,3,4}));

  LatticeGaugeField U(grid);
  U = 1.0;          // FIXME: I am setting to cold configuration
  LatticeColourMatrix g(U.Grid());
  g = 1.0;
  GF_heatbath(U, g, hmc_para.hb_offset, hmc_para.betaMM, hmc_para.table_path, pRNG); //hb_nsweeps before calculate equilibrium value


  // LatticeComplex g_trace(grid);
  // g_trace = trace(g);
  // std::cout << g_trace << std::endl;

  typename LatticeColourMatrix::vector_object::scalar_object g0;
  peekSite(g0, g, Coordinate({0,0,0,0}));
  LatticeColourMatrix g0Inv_x_g = adj(g0) * g;
  LatticeComplex g0Inv_x_g_trace(grid);
  g0Inv_x_g_trace = trace(g0Inv_x_g);
  std::cout << g0Inv_x_g_trace << std::endl;

  std::cout << "Finished!"  << std::endl;
  Grid_finalize();

} // main
