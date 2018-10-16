namespace Grid {
namespace QCD {

// redefined runner function; change HybridMonteCarlo to GF_HybridMonteCarlo
template <class Implementation,
          template <typename, typename, typename> class Integrator,
          class RepresentationsPolicy = NoHirep, class ReaderClass = XmlReader>
class GF_HMCWrapperTemplate: public HMCWrapperTemplate<Implementation, Integrator, RepresentationsPolicy, ReaderClass> {
 public:
  INHERIT_FIELD_TYPES(Implementation);
  typedef Implementation ImplPolicy;  // visible from outside
  template <typename S = NoSmearing<Implementation> >
  using IntegratorType = Integrator<Implementation, S, RepresentationsPolicy>;


  void Run(const HMC_PARA &HMC_para){
    NoSmearing<Implementation> S;
    Runner(S, HMC_para);
  }

  template <class SmearingPolicy>
  void Runner(SmearingPolicy &Smearing, const HMC_PARA &HMC_para) {
    auto UGrid = this->Resources.GetCartesian();
    this->Resources.AddRNGs();
    Field U(UGrid);

    // Can move this outside?
    typedef IntegratorType<SmearingPolicy> TheIntegrator;
    TheIntegrator MDynamics(UGrid, this->Parameters.MD, this->TheAction, Smearing);

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
    if(HMC_para.UInitEquil){
      FieldMetaData header;
      std::string file(HMC_para.UFile);
      // std::string file("./U_softly_fixed_4_M0.5");
      //std::string file("./U_equilibrium");
      NerscIO::readConfiguration(U, header,file);
    }

    //SteepestDescentGaugeFix
    if(HMC_para.SDGF){
      Real alpha=0.1;
      // must be false. cannot be true
      std::cout << "before SDGF dOmegaSquare: " << dOmegaSquare2_no_g(U) << std::endl;
	  LatticeColourMatrix g(U._grid);
  	  g = 1.0;
	  GF_heatbath(U, g, 50, HMC_para.betaMM, HMC_para.hb_multi_hit);
	  SU3::GaugeTransform(U, g);
      //FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(U,alpha,100,1.0e-12, 1.0e-12,false);
      std::cout << "after SDGF dOmegaSquare: " << dOmegaSquare2_no_g(U) << std::endl;
    }

    GF_HybridMonteCarlo<TheIntegrator> HMC(this->Parameters, MDynamics,
                                        this->Resources.GetSerialRNG(),
                                        this->Resources.GetParallelRNG(),
                                        this->Resources.GetObservables(), U);

    // Run it
    std::cout << "-------GF evolve--------" << std::endl;
    HMC.evolve(HMC_para);
  }
};

// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using GF_GenericHMCRunner = GF_HMCWrapperTemplate<PeriodicGimplR, Integrator>;

}  // namespace QCD
}  // namespace Grid
