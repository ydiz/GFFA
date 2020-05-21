namespace Grid{
namespace QCD{

int Nexp = 12;


//implement GF_refresh, S, and update_U
//A inherits from B, B inherits from C; then the constructor in B must be written explicitly.
template <class FieldImplementation, class SmearingPolicy, class RepresentationPolicy>
class GFIntegrator : public Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy> {

public:
INHERIT_FIELD_TYPES(FieldImplementation);

GFIntegrator(GridBase* grid, IntegratorParameters Par,
         ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
    : Integrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
          grid, Par, Aset, Sm) {};

  // RealD GF_S(Field& U, const Momenta_k &KK) {
  //   RealD H;
  //   if(KK.newHp) H = Hp(this->P, KK);
  //   else H = - FieldImplementation::FieldSquareNorm(this->P); //minus sign comes from the fact that P is actually timesI(P).
  //
  //   RealD Hterm;
  //   std::cout << GridLogMessage << "Momentum action H_p = "<<  std::setprecision(8)  << H << "\n";
  //
  //   // Actions
  //   for (int level = 0; level < this->as.size(); ++level) {
  //     for (int actionID = 0; actionID < this->as[level].actions.size(); ++actionID) {
  //       // get gauge field from the SmearingPolicy and
  //       // based on the boolean is_smeared in actionID
  //       Field& Us =
  //           this->Smearer.get_U(this->as[level].actions.at(actionID)->is_smeared);
  //       Hterm = this->as[level].actions.at(actionID)->S(Us);
  //       std::cout << GridLogMessage << "S Level " << level << " term "
  //                 << actionID << " H = " << Hterm << std::endl;
  //       H += Hterm;
  //     }
  //     this->as[level].apply(this->S_hireps, this->Representations, level, H);
  //   }
  //
  //   return H;
  // }

inline void GF_refresh(Field& U, GridParallelRNG& pRNG, const Momenta_k &KK, const HMC_PARA &HMC_para) {
  assert(this->P._grid == U._grid);
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

  // set zero mode to zero // FIXME: is this right?
  std::cout << "Setting zero mode of initial P to zero" << std::endl;
  set_zero_mode_to_zero(this->P);

  this->Smearer.set_Field(U);
  this->Representations.update(U);

  for (int level = 0; level < this->as.size(); ++level) {
    for (int actionID = 0; actionID < this->as[level].actions.size(); ++actionID) {
      Field& Us =
          this->Smearer.get_U(this->as[level].actions.at(actionID)->is_smeared);
      this->as[level].actions.at(actionID)->refresh(Us, pRNG);
    }

    this->as[level].apply(this->refresh_hireps, this->Representations, pRNG);
  }
}

void update_U(Field& U, double ep, const Momenta_k &KK) {
  update_U(this->P, U, ep, KK);

  this->t_U += ep;
  int fl = this->levels - 1;
  std::cout << GridLogIntegrator << "   " << "[" << fl << "] U " << " dt " << ep << " : t_U " << this->t_U << std::endl;
}

void update_U(LatticeGaugeField& Mom, LatticeGaugeField& U, double ep, const Momenta_k &KK) {
  LatticeGaugeField deltaU(Mom._grid);

  if(KK.newHp) deltaU = dHdP(Mom, KK);
  else deltaU = Mom;

  // // // set zero mode to zero // FIXME: Talk to Norman; is doing this right ?
  // set_zero_mode_to_zero(deltaU);
  // // measure_A(deltaU, {{0,0,0,0}, {1,0,0,0}}, false);

  parallel_for(int ss=0;ss<Mom._grid->oSites();ss++){
   for (int mu = 0; mu < Nd; mu++)
     U[ss]._internal[mu] = ProjectOnGroup(Exponentiate(deltaU[ss]._internal[mu], ep, Nexp) * U[ss]._internal[mu]);
  }

  // std::cout << "explicitly setting A(k=0) to 0" << std::endl;
  // set_zero_mode_U_to_zero(U); // FIXME: set zero mode to 0
  // // std::cout << "after setting A(k=0) to 0" << std::endl;
  // // measure_A(U, {{0,0,0,0}});

  this->Smearer.set_Field(U);
  this->Representations.update(U);  // void functions if fundamental representation
}

void GF_integrate(Field& U, const Momenta_k &KK, const HMC_PARA &HMC_para) {
  // reset the clocks
  this->t_U = 0;
  for (int level = 0; level < this->as.size(); ++level) {
    this->t_P[level] = 0;
  }

  for (int step = 0; step < this->Params.MDsteps; ++step) {  // MD step
    int first_step = (step == 0);
    int last_step = (step == this->Params.MDsteps - 1);
    this->step(U, 0, first_step, last_step, KK);

    // For measureing A
    if(HMC_para.measure_A) {
      std::cout << "step: " << step << std::endl;
      std::cout << "measure A(k): " << std::endl;

      using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
      static LatticeUeval last_log_evals(U._grid); 
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

virtual void step(Field& U, int level, int first, int last, const Momenta_k &KK) = 0;

};



}}
