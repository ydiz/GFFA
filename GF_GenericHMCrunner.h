namespace Grid {
namespace QCD {

template <class Implementation,
          template <typename, typename, typename> class Integrator,
          class RepresentationsPolicy = NoHirep, class ReaderClass = XmlReader>
class GF_HMCWrapperTemplate {
 public:
  INHERIT_FIELD_TYPES(Implementation);
  typedef Implementation ImplPolicy;  // visible from outside
  template <typename S = NoSmearing<Implementation> >
  using IntegratorType = Integrator<Implementation, S, RepresentationsPolicy>;

  HMCparameters Parameters;
  std::string ParameterFile;
  HMCResourceManager<Implementation> Resources;

  // The set of actions (keep here for lower level users, for now)
  MyActionSet<Field, RepresentationsPolicy> TheAction;
  GF_HMCWrapperTemplate(const HMCparameters &Par) {
    Parameters = Par;
  }

  void Run(const GFFAParams &HMC_para){
    NoSmearing<Implementation> S;
    Runner(S, HMC_para);
  }

  template <class SmearingPolicy>
  void Runner(SmearingPolicy &Smearing, const GFFAParams &HMC_para) {
    auto UGrid = this->Resources.GetCartesian();
    this->Resources.AddRNGs();
    Field U(UGrid);

    // Can move this outside?
    typedef IntegratorType<SmearingPolicy> TheIntegrator;
    TheIntegrator MDynamics(UGrid, this->Parameters.MD, this->TheAction, Smearing, HMC_para.cell_size);

    if (this->Parameters.StartingType == "HotStart") {
      // Hot start
      this->Resources.SeedFixedIntegers();
      Implementation::HotConfiguration(this->Resources.GetParallelRNG(), U);
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

    Smearing.set_Field(U);

    //read U in equilibrium
    if( !(this->Parameters.StartingType == "CheckpointStart") && !HMC_para.UFile.empty() ){

      this->Resources.SeedFixedIntegers(); // this is enough!
      FieldMetaData header;
      std::string file(HMC_para.UFile);
      // std::string file("./U_softly_fixed_4_M0.5");
      NerscIO::readConfiguration(U, header,file);
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

template <template <typename, typename, typename> class Integrator>
using GF_GenericHMCRunner = GF_HMCWrapperTemplate<PeriodicGimplR, Integrator, NoHirep, JSONReader>;
}  // namespace QCD
}  // namespace Grid
