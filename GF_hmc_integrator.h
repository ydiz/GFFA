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

  RealD GF_S(Field& U, const Momenta_k &KK) {
    RealD H;
    if(KK.newHp) H = Hp(this->P, KK);
    else H = - FieldImplementation::FieldSquareNorm(this->P); //minus sign comes from the fact that P is actually timesI(P).

    RealD Hterm;
    std::cout << GridLogMessage << "Momentum action H_p = " << H << "\n";

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
      this->as[level].apply(this->S_hireps, this->Representations, level, H);
    }

    return H;
  }



inline void GF_refresh(Field& U, GridParallelRNG& pRNG, const Momenta_k &KK) {
  assert(this->P._grid == U._grid);
  std::cout << GridLogIntegrator << "Integrator refresh\n";

  if(KK.newHp) GF_generate_P(this->P, pRNG, KK);
  else FieldImplementation::generate_momenta(this->P, pRNG);

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

  parallel_for(int ss=0;ss<Mom._grid->oSites();ss++){
   for (int mu = 0; mu < Nd; mu++)
     U[ss]._internal[mu] = ProjectOnGroup(Exponentiate(deltaU[ss]._internal[mu], ep, Nexp) * U[ss]._internal[mu]);
  }

  this->Smearer.set_Field(U);
  this->Representations.update(U);  // void functions if fundamental representation
}

void GF_integrate(Field& U, const Momenta_k &KK) {
  // reset the clocks
  this->t_U = 0;
  for (int level = 0; level < this->as.size(); ++level) {
    this->t_P[level] = 0;
  }

  for (int step = 0; step < this->Params.MDsteps; ++step) {  // MD step
    int first_step = (step == 0);
    int last_step = (step == this->Params.MDsteps - 1);
    this->step(U, 0, first_step, last_step, KK);
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



template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class GFLeapFrog : public GFIntegrator<FieldImplementation, SmearingPolicy,
                                   RepresentationPolicy> {
 public:
  typedef GFLeapFrog<FieldImplementation, SmearingPolicy, RepresentationPolicy>
      Algorithm;
  INHERIT_FIELD_TYPES(FieldImplementation);

  std::string integrator_name(){return "GFLeapFrog";}

  GFLeapFrog(GridBase* grid, IntegratorParameters Par,
           ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
      : GFIntegrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm){};

  void step(Field& U, int level, int first, int last) {}

  void step(Field& U, int level, int _first, int _last, const Momenta_k &KK)
  {
    int fl = this->as.size() - 1;

    RealD eps = this->Params.trajL/this->Params.MDsteps;
    for (int l = 0; l <= level; ++l) eps /= this->as[l].multiplier;

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {
      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, eps / 2.0);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, last_step, KK);
      }

      int mm = last_step ? 1 : 2;
      this->update_P(U, level, mm * eps / 2.0);

    }
  }
};


}}
