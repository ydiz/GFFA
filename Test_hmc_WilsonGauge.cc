#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
  using namespace Grid::QCD;

  Grid_init(&argc, &argv);
  GridLogLayout();

   // Typedefs to simplify notation
  typedef GenericHMCRunner<LeapFrog> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  // typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  // TopologyObsParameters TopParams;
  // TopParams.interval = 5;
  // TopParams.do_smearing = true;
  // TopParams.Smearing.steps = 200;
  // TopParams.Smearing.step_size = 0.01;
  // TopParams.Smearing.meas_interval = 50;
  // TopParams.Smearing.maxTau = 2.0;
  // TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  RealD beta = 5.6;
  WilsonGaugeActionR Waction(beta);

  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  //Level1.push_back(WGMod.getPtr());
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
