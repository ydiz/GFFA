#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef WilsonImplR FermionImplPolicy;
  typedef MobiusFermionR FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IntegratorParameters MD;
  typedef GenericHMCRunner<LeapFrog> HMCWrapper; 
  MD.name    = std::string("Leap Frog");
  //  typedef GenericHMCRunner<ForceGradient> HMCWrapper; 
  //  MD.name    = std::string("Force Gradient");
  //typedef GenericHMCRunner<MinimumNorm2> HMCWrapper; 
  //  MD.name    = std::string("MinimumNorm2");
  MD.MDsteps = 40;
  MD.trajL   = 1.0;
  
  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 0;
  HMCparams.Trajectories     = 1;
  HMCparams.NoMetropolisUntil=  20;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.MD = MD;
  HMCWrapper TheHMC(HMCparams);

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_EODWF_lat";
  CPparams.rng_prefix    = "ckpoint_EODWF_rng";
  CPparams.saveInterval  = 10;
  CPparams.format        = "IEEE64BIG";
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  const int Ls      = 4;
  Real beta         = 2.13;
  Real light_mass   = 0.01;
  // Real strange_mass = 0.04;
  Real pv_mass      = 1.0;
  RealD M5  = 1.8;
  RealD b   = 1.5; // Scale factor two
  RealD c   = 0.5;

  // RHMC
  // OneFlavourRationalParams OFRp;
  // OFRp.lo       = 1.0e-2;
  // OFRp.hi       = 25;
  // OFRp.MaxIter  = 10000;
  // OFRp.tolerance= 1.0e-7;
  // OFRp.degree   = 10;
  // OFRp.precision= 40;

  // EOFA
  OneFlavourRationalParams OFRp;
  OFRp.lo       = 0.98;
  OFRp.hi       = 25.0;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-7;
  OFRp.degree   = 10;
  OFRp.precision= 40;

  std::vector<Real> hasenbusch({ 0.1 });

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  IwasakiGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  LatticeGaugeField U(GridPtr);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams Params(boundary);
  
  double StoppingCondition = 1e-10;
  double MaxCGIterations = 30000;
  ConjugateGradient<LatticeFermion>  CG(StoppingCondition,MaxCGIterations);

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field> Level1(1);
  ActionLevel<HMCWrapper::Field> Level2(4);

  ////////////////////////////////////
  // Strange action
  ////////////////////////////////////

  // Setup for RHMC
  //  FermionAction StrangeOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_mass,M5,b,c, Params);
  //  OneFlavourRationalPseudoFermionAction<FermionImplPolicy> StrangePseudoFermion(StrangeOp,OFRp);
  // Level1.push_back(&StrangePseudoFermion);
  
  // DJM: setup for EOFA ratio (Shamir)
  DomainWallEOFAFermionR Strange_Op_L(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, strange_mass, strange_mass, pv_mass, 0.0, -1, M5);
  DomainWallEOFAFermionR Strange_Op_R(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, pv_mass, strange_mass, pv_mass, -1.0, 1, M5);
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> EOFA(Strange_Op_L, Strange_Op_R, CG, OFRp, true);
  Level1.push_back(&EOFA);

  // // DJM: setup for EOFA ratio (Mobius)
  // MobiusEOFAFermionR Strange_Op_L(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, strange_mass, strange_mass, pv_mass, 0.0, -1, M5, b, c);
  // MobiusEOFAFermionR Strange_Op_R(U, *FGrid, *FrbGrid, *GridPtr, *GridRBPtr, pv_mass, strange_mass, pv_mass, -1.0, 1, M5, b, c);
  // ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy> EOFA(Strange_Op_L, Strange_Op_R, CG, OFRp, true);
  // Level1.push_back(&EOFA);

  // Setup for RHMC ratio
  // FermionAction StrangeOp (U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,strange_mass,M5,b,c, Params);
  // FermionAction StrangePauliVillarsOp(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,pv_mass,  M5,b,c, Params);
  // OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> StrangePseudoFermion(StrangePauliVillarsOp,StrangeOp,OFRp);
  // Level1.push_back(&StrangePseudoFermion);

  ////////////////////////////////////
  // up down action
  ////////////////////////////////////
  std::vector<Real> light_den;
  std::vector<Real> light_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(light_mass);
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
    light_num.push_back(hasenbusch[h]);
  }
  light_num.push_back(pv_mass);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;

  for(int h=0;h<n_hasenbusch+1;h++){
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;
    Numerators.push_back  (new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_num[h],M5,b,c, Params));
    Denominators.push_back(new FermionAction(U,*FGrid,*FrbGrid,*GridPtr,*GridRBPtr,light_den[h],M5,b,c, Params));
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],CG,CG));
  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Level1.push_back(Quotients[h]);
  }

  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main
