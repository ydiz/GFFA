namespace Grid {

template <class Implementation,
          template <typename, typename> class Integrator,
          class RepresentationsPolicy = NoHirep, class ReaderClass = XmlReader>
class GF_HMCWrapperTemplate {
 public:
  INHERIT_FIELD_TYPES(Implementation);
  typedef Implementation ImplPolicy;  // visible from outside
  using TheIntegrator = Integrator<Implementation, RepresentationsPolicy>;

  HMCparameters Parameters;
  std::string ParameterFile;
  HMCResourceManager<Implementation> Resources;

  // The set of actions (keep here for lower level users, for now)
  MyActionSet<Field, RepresentationsPolicy> TheAction;
  GF_HMCWrapperTemplate(const HMCparameters &Par) {
    Parameters = Par;
  }

  void Run(const GFFAParams &HMC_para){
  //   NoSmearing<Implementation> S;
  //   Runner(S, HMC_para);
  // }
  //
  // template <class SmearingPolicy>
  // void Runner(SmearingPolicy &Smearing, const GFFAParams &HMC_para) {
    auto UGrid = this->Resources.GetCartesian();
    this->Resources.AddRNGs();
    Field U(UGrid);

    TheIntegrator MDynamics(UGrid, this->Parameters.MD, this->TheAction, HMC_para.cell_size);

    if (this->Parameters.StartingType == "HotStart") {
      // Hot start
      this->Resources.SeedFixedIntegers();
      // Implementation::HotConfiguration(this->Resources.GetParallelRNG(), U); // the scale of LieRandomize is 1.0 // Grid HotStart is not really hot start

      LatticeColourMatrix Umu(U.Grid());
      for (int mu = 0; mu < Nd; mu++) {
        SU3::LieRandomize(this->Resources.GetParallelRNG(), Umu, 2. * M_PI);
        PokeIndex<LorentzIndex>(U, Umu, mu);
      }

      std::cout << "Hot Starting Configuration: Plaq " << WilsonLoops<PeriodicGimplR>::avgPlaquette(U) << std::endl;
      std::cout << "Hot Starting Configuration: LinkTrace " << WilsonLoops<PeriodicGimplR>::linkTrace(U) << std::endl;

    } else if (this->Parameters.StartingType == "ColdStart") {
      // Cold start
      this->Resources.SeedFixedIntegers();
      Implementation::ColdConfiguration(this->Resources.GetParallelRNG(), U);
    } else if (this->Parameters.StartingType == "TepidStart") {
      // Tepid start
      this->Resources.SeedFixedIntegers();
      Implementation::TepidConfiguration(this->Resources.GetParallelRNG(), U);
    } else if (this->Parameters.StartingType == "CheckpointStart") {
      // CheckpointRestart
      this->Resources.GetCheckPointer()->CheckpointRestore(this->Parameters.StartTrajectory, U,
                                   this->Resources.GetSerialRNG(),
                                   this->Resources.GetParallelRNG());
    }

    //read U in equilibrium
    if( !(this->Parameters.StartingType == "CheckpointStart") && !HMC_para.UFile.empty() ){

      this->Resources.SeedFixedIntegers(); // this is enough!
      FieldMetaData header;
      NerscIO::readConfiguration(U, header, HMC_para.UFile);
    }

    GF_HybridMonteCarlo<TheIntegrator> HMC(this->Parameters, MDynamics,
                                        this->Resources.GetSerialRNG(),
                                        this->Resources.GetParallelRNG(),
                                        this->Resources.GetObservables(), U);

    // Run it
    std::cout << "-------GF evolve--------" << std::endl;
    std::cout << "Integrator name: " << MDynamics.integrator_name() << std::endl;
    HMC.evolve(HMC_para);
  }
};

// template <template <typename, typename, typename> class Integrator>
// using GF_GenericHMCRunner = GF_HMCWrapperTemplate<PeriodicGimplR, Integrator>;

// template <template <typename, typename, typename> class Integrator>
// using GF_GenericHMCRunner = GF_HMCWrapperTemplate<PeriodicGimplR, Integrator, NoHirep, JSONReader>;
//
template <template <typename, typename> class Integrator>
using GF_GenericHMCRunner = GF_HMCWrapperTemplate<PeriodicGimplR, Integrator, NoHirep, JSONReader>;

}  // namespace Grid
