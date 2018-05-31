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
               ActionSet<Field, RepresentationPolicy>& Aset, SmearingPolicy& Sm)
      : GFIntegrator<FieldImplementation, SmearingPolicy, RepresentationPolicy>(
            grid, Par, Aset, Sm){};

  std::string integrator_name(){return "GFMininumNorm2";}

  void step(Field& U, int level, int first, int last) {}

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










}
}
