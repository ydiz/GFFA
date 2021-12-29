#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonImplR FermionImplPolicy;
  typedef WilsonFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 5;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  // RealD beta = 5.6 ;  // ?? which one is more common
  Real beta         = 2.13;
  int Ls = 4;
  Real mass = 0.3;
  Real M5 = 1.8;

    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);
  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);


  // Can we define an overloaded operator that does not need U and initialises
  // it with zeroes?
  // FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);
  DomainWallFermionR fermOp(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, mass, M5);
  ConjugateGradient<FermionField> CG(1.0e-8, 2000);
  TwoFlavourPseudoFermionAction<FermionImplPolicy> Nf2(fermOp, CG, CG);
  Nf2.is_smeared = false;


    // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Nf2);
  TheHMC.TheAction.push_back(Level1);

  WilsonGaugeActionR Waction(beta);
  ActionLevel<HMCWrapper::Field> Level2(4);
  Level2.push_back(&Waction);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////


  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main
