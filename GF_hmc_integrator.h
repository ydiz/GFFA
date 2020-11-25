namespace Grid{
namespace QCD{

int Nexp = 12;


//implement GF_refresh, S, and update_U
//A inherits from B, B inherits from C; then the constructor in B must be written explicitly.
template <class FieldImplementation, class SmearingPolicy, class RepresentationPolicy>
// class GFIntegrator : public Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy> {
class GFIntegrator  {

public:
// INHERIT_FIELD_TYPES(FieldImplementation);

// typedef typename FieldImplementation::Field MomentaField;  //for readability

// GFIntegrator(GridBase* grid, IntegratorParameters Par,
//          ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
//     : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
//           grid, Par, Aset, Sm) {};

// const MyActionSet<Field, RepresentationPolicy> as; // replace the ActionSet in Integrator class

  typedef typename FieldImplementation::Field MomentaField;  //for readability
  typedef typename FieldImplementation::Field Field;

  int levels;  // number of integration levels
  double t_U;  // Track time passing on each level and for U and for P
  std::vector<double> t_P;  

  MomentaField P;
  SmearingPolicy& Smearer;
  RepresentationPolicy Representations;
  IntegratorParameters Params;

  const MyActionSet<Field, RepresentationPolicy> as;


GFIntegrator(GridBase* grid, IntegratorParameters Par,
           MyActionSet<Field, RepresentationPolicy>& Aset,
           SmearingPolicy& Sm)
  : Params(Par),
    as(Aset),
    P(grid),
    levels(Aset.size()),
    Smearer(Sm),
    Representations(grid) 
{
  t_P.resize(levels, 0.0);
  t_U = 0.0;
  // initialization of smearer delegated outside of Integrator
};



RealD GF_S(Field& U, const Momenta_k &KK) {
  RealD H;
  if(KK.newHp) H = Hp(this->P, KK);
  else H = - FieldImplementation::FieldSquareNorm(this->P) / HMC_MOMENTUM_DENOMINATOR; // - trace (P*P) / denom // minus sign comes from the fact that P is actually timesI(P).  // ! Must add HMC_MOMENTUM_DENOMINATOR factor for new version Grid

  RealD Hterm;
  std::cout << GridLogMessage << "Momentum action H_p = "<<  std::setprecision(8)  << H << "\n";

  // Actions
  for (int level = 0; level < this->as.size(); ++level) {
    for (int actionID = 0; actionID < this->as[level].actions.size(); ++actionID) {
      // get gauge field from the SmearingPolicy and
      // based on the boolean is_smeared in actionID
      Field& Us =
          this->Smearer.get_U(this->as[level].actions.at(actionID)->is_smeared);
      Hterm = this->as[level].actions.at(actionID)->S(Us);
      std::cout << GridLogMessage << "S Level " << level << " term "
                << actionID << " H = " << Hterm << std::endl;
      H += Hterm;
    }
    // this->as[level].apply(this->S_hireps, this->Representations, level, H);
  }

  return H;
}

inline void GF_refresh(Field& U, GridParallelRNG& pRNG, const Momenta_k &KK, const GFFAParams &HMC_para) {
  assert(this->P.Grid() == U.Grid());
  std::cout << GridLogIntegrator << "Integrator refresh\n";

  // if(HMC_para.measure_A) {
  //   double fixed_P_k = HMC_para.fixed_P_k; // start with P(k)^a = fixed_P_k expect for P(k=0)^a = 0
  //   //
  //   // double P_n0 = fixed_P_k * std::sqrt(KK.vol);
  //   //
  //   // this->P = 0.;
  //   // LatticeGaugeField::vector_object::scalar_object P_site0;
  //   //
  //   // SU3::Matrix tmp; tmp = zero;
  //   // for(int a=0; a<8; ++a) {
  //   //   SU3::Matrix ta;
  //   //   SU3::generator(a, ta);
  //   //   tmp = tmp + ta * P_n0;
  //   // }
  //   // for (int mu = 0; mu < 4; mu++) P_site0(mu) = tmp();
  //   // pokeSite(timesI(P_site0), this->P, {0,0,0,0});
  //
  //   // All site of Pk should be the same
  //   std::cout << "Initial momenta: P_mu(k)^a = " + std::to_string(fixed_P_k) + ", not random" << std::endl;
  //
  //   FieldImplementation::generate_momenta(this->P, pRNG);
  // }
  // else {
  if(KK.newHp) GF_generate_P(this->P, pRNG, KK);
  else FieldImplementation::generate_momenta(this->P, pRNG);
  // }

  // // I am not setting initial zero mode P(k=0) to zero
  // // set zero mode to zero // FIXME: is this right?
  // if(KK.newHp) {
  //   std::cout << "Setting zero mode of initial P to zero" << std::endl;
  //   set_zero_mode_to_zero(this->P);
  // }

  this->Smearer.set_Field(U);
  this->Representations.update(U);

  for (int level = 0; level < this->as.size(); ++level) {
    for (int actionID = 0; actionID < this->as[level].actions.size(); ++actionID) {
      Field& Us =
          this->Smearer.get_U(this->as[level].actions.at(actionID)->is_smeared);
      this->as[level].actions.at(actionID)->refresh(Us, pRNG);
    }

    // this->as[level].apply(this->refresh_hireps, this->Representations, pRNG);
  }
}

void update_U(Field& U, double ep, const Momenta_k &KK) {
  update_U(this->P, U, ep, KK);

  this->t_U += ep;
  int fl = this->levels - 1;
  std::cout << GridLogIntegrator << "   " << "[" << fl << "] U " << " dt " << ep << " : t_U " << this->t_U << std::endl;
}

void update_U(LatticeGaugeField& Mom, LatticeGaugeField& U, double ep, const Momenta_k &KK) {
  LatticeGaugeField deltaU(Mom.Grid());

  if(KK.newHp) deltaU = dHdP(Mom, KK);
  else deltaU = Mom;

  // // // set zero mode to zero // FIXME: Talk to Norman; is doing this right ?
  // set_zero_mode_to_zero(deltaU);
  // // measure_A(deltaU, {{0,0,0,0}, {1,0,0,0}}, false);

  // auto U_v = U.View();
  // auto deltaU_v = deltaU.View();
  autoView(U_v, U, AcceleratorWrite);
  autoView(deltaU_v, deltaU, AcceleratorRead);

  // parallel_for(int ss=0;ss<Mom.Grid()->oSites();ss++){
  // thread_for(ss, Mom.Grid()->oSites(), {
  accelerator_for(ss, Mom.Grid()->oSites(), vComplex::Nsimd(), {
   for (int mu = 0; mu < Nd; mu++)
     // U[ss]._internal[mu] = ProjectOnGroup(Exponentiate(deltaU[ss]._internal[mu], ep, Nexp) * U[ss]._internal[mu]);
     U_v[ss](mu) = ProjectOnGroup(Exponentiate(deltaU_v[ss](mu), ep, Nexp) * U_v[ss](mu));
  });

  // std::cout << "explicitly setting A(k=0) to 0" << std::endl;
  // set_zero_mode_U_to_zero(U); // FIXME: set zero mode to 0
  // // std::cout << "after setting A(k=0) to 0" << std::endl;
  // // measure_A(U, {{0,0,0,0}});

  this->Smearer.set_Field(U);
  this->Representations.update(U);  // void functions if fundamental representation
}

void update_P(Field& U, int level, double ep, GridSerialRNG &sRNG, GridParallelRNG &pRNG) 
{
  this->t_P[level] += ep;
  update_P(this->P, U, level, ep, sRNG, pRNG);

  std::cout << GridLogIntegrator << "[" << level << "] P " << " dt " << ep << " : t_P " << this->t_P[level] << std::endl;
}



void update_P(MomentaField& Mom, Field& U, int level, double ep, GridSerialRNG &sRNG, GridParallelRNG &pRNG) {
  // input U actually not used in the fundamental case
  // Fundamental updates, include smearing

  for (int a = 0; a < this->as[level].actions.size(); ++a) {
    double start_full = usecond();
    Field force(U.Grid());
    conformable(U.Grid(), Mom.Grid());

    Field& Us = this->Smearer.get_U(this->as[level].actions.at(a)->is_smeared);
    double start_force = usecond();
    this->as[level].actions.at(a)->deriv(Us, force, sRNG, pRNG);  // deriv should NOT include Ta // zyd: should still be fine if deriv has Ta when smearing is off

    std::cout << GridLogIntegrator << "Smearing (on/off): " << this->as[level].actions.at(a)->is_smeared << std::endl;
    if (this->as[level].actions.at(a)->is_smeared) this->Smearer.smeared_force(force);
    force = FieldImplementation::projectForce(force); // Ta for gauge fields
    double end_force = usecond();
    Real force_abs = std::sqrt(norm2(force)/U.Grid()->gSites());
    std::cout << GridLogIntegrator << "["<<level<<"]["<<a<<"] Force average: " << force_abs << std::endl;
    Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
    double end_full = usecond();
    double time_full  = (end_full - start_full) / 1e3;
    double time_force = (end_force - start_force) / 1e3;
    std::cout << GridLogMessage << "["<<level<<"]["<<a<<"] P update elapsed time: " << time_full << " ms (force: " << time_force << " ms)"  << std::endl;
  }

  // // Force from the other representations
  // as[level].apply(update_P_hireps, Representations, Mom, U, ep);
}





void GF_integrate(Field& U, const Momenta_k &KK, const GFFAParams &HMC_para, GridSerialRNG &sRNG, GridParallelRNG &pRNG) {
  // reset the clocks
  this->t_U = 0;
  for (int level = 0; level < this->as.size(); ++level) {
    this->t_P[level] = 0;
  }

  for (int step = 0; step < this->Params.MDsteps; ++step) {  // MD step
    int first_step = (step == 0);
    int last_step = (step == this->Params.MDsteps - 1);
    this->step(U, 0, first_step, last_step, KK, sRNG, pRNG);

    // For measureing A
    if(HMC_para.measure_A) {
      std::cout << "step: " << step << std::endl;
      std::cout << "measure A(k): " << std::endl;

      using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
      static LatticeUeval last_log_evals(U.Grid()); 
      static bool last_log_evals_initialized = false;

      measure_A(U, HMC_para.measure_A_coors, last_log_evals, last_log_evals_initialized);
    }
  }

  // Check the clocks all match on all levels
  for (int level = 0; level < this->as.size(); ++level) {
    assert(fabs(this->t_U - this->t_P[level]) < 1.0e-6);  // must be the same
    std::cout << GridLogIntegrator << " times[" << level
              << "]= " << this->t_P[level] << " " << this->t_U << std::endl;
  }

  // and that we indeed got to the end of the trajectory
  assert(fabs(this->t_U - this->Params.trajL) < 1.0e-6);

}

virtual void step(Field& U, int level, int first, int last, const Momenta_k &KK, GridSerialRNG &sRNG, GridParallelRNG &pRNG) = 0;




  virtual ~GFIntegrator() {}

  virtual std::string integrator_name() = 0;

  void print_parameters()
  {
    std::cout << GridLogMessage << "[Integrator] Name : "<< integrator_name() << std::endl;
    Params.print_parameters();
  }

  void print_actions()
  {
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
    std::cout << GridLogMessage << "[Integrator] Action summary: "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      std::cout << GridLogMessage << "[Integrator] ---- Level: "<< level << std::endl;
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
  std::cout << GridLogMessage << "["<< as[level].actions.at(actionID)->action_name() << "] ID: " << actionID << std::endl;
  std::cout << as[level].actions.at(actionID)->LogParameters();
      }
    }
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::"<< std::endl;

  }




};



}}
