namespace Grid{
namespace QCD{

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
           MyActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
      : GFIntegrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm){};

  void step(Field& U, int level, int first, int last) {assert(0);}

  void step(Field& U, int level, int _first, int _last, const Momenta_k &KK, GridSerialRNG &sRNG, GridParallelRNG &pRNG)
  {
    int fl = this->as.size() - 1;

    RealD eps = this->Params.trajL/this->Params.MDsteps;
    for (int l = 0; l <= level; ++l) eps /= this->as[l].multiplier;

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {
      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        // this->update_P(U, level, eps / 2.0, sRNG, pRNG);
        this->update_P(U, level, eps / 2.0, sRNG, pRNG, true);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, last_step, KK, sRNG, pRNG);
      }

      int mm = last_step ? 1 : 2;
      // this->update_P(U, level, mm * eps / 2.0, sRNG, pRNG);
      this->update_P(U, level, mm * eps / 2.0, sRNG, pRNG, false);

    }
  }
};


template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class GFMinimumNorm2 : public GFIntegrator<FieldImplementation, SmearingPolicy,
                                       RepresentationPolicy> {
 private:
  const RealD lambda = 0.1931833275037836;

 public:
  INHERIT_FIELD_TYPES(FieldImplementation);

  GFMinimumNorm2(GridBase* grid, IntegratorParameters Par,
               MyActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
      : GFIntegrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm){};

  std::string integrator_name(){return "GFMininumNorm2";}

  void step(Field& U, int level, int first, int last) {assert(0);}

  void step(Field& U, int level, int _first, int _last, const Momenta_k &KK) {

    int fl = this->as.size() - 1;

    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    int multiplier = this->as[level].multiplier;
    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, 0, KK);
      }

      this->update_P(U, level, (1.0 - 2.0 * lambda) * eps);

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, 0, last_step, KK);
      }

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }
  }
};

// add KK to update_U() and step()
template <class FieldImplementation, class SmearingPolicy,
          class RepresentationPolicy =
              Representations<FundamentalRepresentation> >
class GFForceGradient : public GFIntegrator<FieldImplementation, SmearingPolicy,
                                       RepresentationPolicy> {
  private:
  const RealD lambda = 1.0 / 6.0;
  const RealD chi = 1.0 / 72.0;
  const RealD xi = 0.0;
  const RealD theta = 0.0;

 public:
  INHERIT_FIELD_TYPES(FieldImplementation);

  GFForceGradient(GridBase* grid, IntegratorParameters Par,
               MyActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
      : GFIntegrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm){};

  std::string integrator_name(){return "GFForceGradient";}

  void FG_update_P(Field& U, int level, double fg_dt, double ep, const Momenta_k &KK) {
    Field Ufg(U.Grid());
    Field Pfg(U.Grid());
    Ufg = U;
    Pfg = Zero();
    std::cout << GridLogIntegrator << "FG update " << fg_dt << " " << ep
              << std::endl;
    // prepare_fg; no prediction/result cache for now
    // could relax CG stopping conditions for the
    // derivatives in the small step since the force gets multiplied by
    // a tiny dt^2 term relative to main force.
    //
    // Presently 4 force evals, and should have 3, so 1.33x too expensive.
    // could reduce this with sloppy CG to perhaps 1.15x too expensive
    // even without prediction.
    this->update_P(Pfg, Ufg, level, 1.0);
    this->update_U(Pfg, Ufg, fg_dt, KK);
    this->update_P(Ufg, level, ep);
  }

  void step(Field& U, int level, int first, int last) { assert(0);}

  void step(Field& U, int level, int _first, int _last, const Momenta_k &KK) {

    RealD eps = this->Params.trajL/this->Params.MDsteps * 2.0;
    for (int l = 0; l <= level; ++l) eps /= 2.0 * this->as[l].multiplier;

    RealD Chi = chi * eps * eps * eps;

    int fl = this->as.size() - 1;

    int multiplier = this->as[level].multiplier;

    for (int e = 0; e < multiplier; ++e) {  // steps per step

      int first_step = _first && (e == 0);
      int last_step = _last && (e == multiplier - 1);

      if (first_step) {  // initial half step
        this->update_P(U, level, lambda * eps);
      }

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, first_step, 0, KK);
      }

      this->FG_update_P(U, level, 2 * Chi / ((1.0 - 2.0 * lambda) * eps),
                        (1.0 - 2.0 * lambda) * eps, KK);

      if (level == fl) {  // lowest level
        this->update_U(U, 0.5 * eps, KK);
      } else {  // recursive function call
        this->step(U, level + 1, 0, last_step, KK);
      }

      int mm = (last_step) ? 1 : 2;
      this->update_P(U, level, lambda * eps * mm);
    }

  }

};







}
}
