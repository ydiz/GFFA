namespace Grid{

// int Nexp = 12;
int Nexp = 20; // Taylor expansion with Nexp=20 is very close to Cayley-Hamilton


//implement GF_refresh, S, and update_U
//A inherits from B, B inherits from C; then the constructor in B must be written explicitly.
template <class FieldImplementation, class RepresentationPolicy>
class GFIntegrator  {

public:
  typedef typename FieldImplementation::Field MomentaField;  //for readability
  typedef typename FieldImplementation::Field Field;

  int levels;  // number of integration levels
  double t_U;  // Track time passing on each level and for U and for P
  std::vector<double> t_P;  

  MomentaField P;
  IntegratorParameters Params;

  const MyActionSet<Field, RepresentationPolicy> as;

  LatticeLorentzScalar mask; // For update only a cell

GFIntegrator(GridBase* grid, IntegratorParameters Par,
           MyActionSet<Field, RepresentationPolicy>& Aset,
           const Coordinate &cell_size)
  : Params(Par),
    as(Aset),
    P(grid),
    levels(Aset.size()),
    mask(grid) // For update only a cell
{
  t_P.resize(levels, 0.0);
  t_U = 0.0;

  mask = get_mask(grid, cell_size);
};



RealD GF_S(Field& U, const Momenta_k &KK) {
  RealD H;
  if(KK.newHp) H = Hp(this->P, KK);
  else H = - FieldImplementation::FieldSquareNorm(this->P) / HMC_MOMENTUM_DENOMINATOR; // - trace (P*P) / denom // minus sign comes from the fact that P is actually timesI(P).  // ! Must add HMC_MOMENTUM_DENOMINATOR factor for new version Grid

  RealD Hterm;
  std::cout << GridLogMessage << "Momentum action H_p = "<<  std::setprecision(8)  << H << "\n";

  for (int level = 0; level < this->as.size(); ++level) {
    for (int actionID = 0; actionID < this->as[level].actions.size(); ++actionID) {

      Hterm = this->as[level].actions.at(actionID)->S(U);
      std::cout << GridLogMessage << "S Level " << level << " term "
                << actionID << " H = " << Hterm << std::endl;
      H += Hterm;
    }
  }

  return H;
}

inline void GF_refresh(Field& U, GridSerialRNG & sRNG, GridParallelRNG& pRNG, const Momenta_k &KK, const GFFAParams &HMC_para) {
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
  if(KK.newHp) {
    if(HMC_para.isCell) {
      // std::cout << "before cell_grid" << std::endl;
      GridBase *cell_grid = pRNG.Grid();
      Coordinate cell_size = cell_grid->_fdimensions;

      LatticeGaugeField tmp(cell_grid);
      GF_generate_P(tmp, pRNG, KK);

      localCopyRegion(tmp, this->P, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
      this->P = this->P * mask; // set "links that are perpendicular to surface" to 0.
      // std::cout << "after cell_grid" << std::endl;
    }
    else GF_generate_P(this->P, pRNG, KK);
  }
  else {
    if(HMC_para.isCell) {
      GridBase *cell_grid = pRNG.Grid();
      Coordinate cell_size = cell_grid->_fdimensions;

      LatticeGaugeField tmp(cell_grid);
      FieldImplementation::generate_momenta(tmp, pRNG);

      localCopyRegion(tmp, this->P, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
      this->P = this->P * mask; // set "links that are perpendicular to surface" to 0.
    }
    else FieldImplementation::generate_momenta(this->P, pRNG);
  }
  // }

  for (int level = 0; level < as.size(); ++level) {
    for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
      // get gauge field from the SmearingPolicy and
      // based on the boolean is_smeared in actionID
      // Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
      as[level].actions.at(actionID)->refresh(U, sRNG, pRNG); // This is necessary when adding fermion action
    }
    // Refresh the higher representation actions
    // as[level].apply(refresh_hireps, Representations, sRNG, pRNG);
  }

  // std::cout << "Not doing gauge transformation at the beginning of a trajectory."  << std::endl;

  static bool first_traj = true;
  if(HMC_para.isGFFA && HMC_para.isCell && first_traj) {  // gauge transformation must be done only before the first trajectory!
    std::cout << "Doing gauge transformation at the beginning of a trajectory."  << std::endl;

    GridBase *cell_grid = KK.Grid();
    Coordinate cell_size = cell_grid->_fdimensions;

    LatticeGaugeField U_cell(cell_grid); U_cell = Zero();
    localCopyRegion(U, U_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
    LatticeLorentzScalar cell_mask = get_cell_mask(cell_grid);
    U_cell = U_cell * cell_mask;
    // Gauge transform U inside the cell
    LatticeColourMatrix g_cell(cell_grid); g_cell = 1.0;
    GridParallelRNG pRNG_cell(cell_grid);      pRNG_cell.SeedFixedIntegers(std::vector<int>({45,12,81,9}));
    GF_heatbath(U_cell, g_cell, HMC_para.innerMC_N, HMC_para.betaMM, HMC_para.table_path, pRNG_cell); //hb_nsweeps before calculate equilibrium value

    LatticeColourMatrix g(U.Grid()); g = 1.0;
    localCopyRegion(g_cell, g, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

    SU<3>::GaugeTransform(U, g);
    std::cout << "Plaquette after gauge transformation: "<< WilsonLoops<PeriodicGimplR>::avgPlaquette(U) << std::endl;
    std::cout << "Link trace after gauge transformation: "<< WilsonLoops<PeriodicGimplR>::linkTrace(U) << std::endl;

    first_traj = false;
  }
  // std::cout << "end of GF_refresh" << std::endl;

}

void update_U(Field& U, double ep, const Momenta_k &KK) {
  update_U(this->P, U, ep, KK);

  this->t_U += ep;
  int fl = this->levels - 1;
}

void update_U(LatticeGaugeField& Mom, LatticeGaugeField& U, double ep, const Momenta_k &KK) {
  // std::cout << "before update_U" << std::endl;
  LatticeGaugeField deltaU(Mom.Grid());

  if(KK.newHp) {
    if(KK.isCell) {
      GridBase *cell_grid = KK.Grid();
      Coordinate cell_size = cell_grid->_fdimensions;

      LatticeGaugeField Mom_cell(cell_grid);
      localCopyRegion(Mom, Mom_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

      LatticeGaugeField tmp(cell_grid);
      tmp = dHdP(Mom_cell, KK);

      localCopyRegion(tmp, deltaU, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
      deltaU = deltaU * mask; // set "links that are perpendicular to surface" to 0.
    }
    else deltaU = dHdP(Mom, KK);
  }
  else {
    deltaU = Mom;
    if(KK.isCell) deltaU = deltaU * mask; // set "links that are perpendicular to surface" to 0.
  }


  autoView(U_v, U, AcceleratorWrite);
  autoView(deltaU_v, deltaU, AcceleratorRead);
  accelerator_for(ss, Mom.Grid()->oSites(), vComplex::Nsimd(), {
   for (int mu = 0; mu < Nd; mu++)
     U_v[ss](mu) = ProjectOnGroup(Exponentiate(deltaU_v[ss](mu), ep, Nexp) * U_v[ss](mu));
  });
  // std::cout << "after update_U" << std::endl;
}

void update_P(Field& U, int level, double ep, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) 
{
  this->t_P[level] += ep;
  update_P(this->P, U, level, ep, sRNG, pRNG, first_step);
}



void update_P(MomentaField& Mom, Field& U, int level, double ep, GridSerialRNG &sRNG, GridParallelRNG &pRNG, bool first_step) {
  // std::cout << "before update_P" << std::endl;
  conformable(U.Grid(), Mom.Grid());

  for (int a = 0; a < this->as[level].actions.size(); ++a) {
    Field force(U.Grid());

    this->as[level].actions.at(a)->deriv(U, force, sRNG, pRNG, first_step);  // deriv should NOT include Ta // zyd: should still be fine if deriv has Ta when smearing is off

    force = FieldImplementation::projectForce(force); // same as Ta(force) 
    Mom -= force * ep * HMC_MOMENTUM_DENOMINATOR; 
  // std::cout << "after update_P" << std::endl;
  }
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
  }
  // and that we indeed got to the end of the trajectory
  assert(fabs(this->t_U - this->Params.trajL) < 1.0e-6);


//   { // calculate average plaq inside the cell
//     GridBase *cell_grid = KK.Grid();
//     Coordinate cell_size = cell_grid->_fdimensions;
//
//     LatticeGaugeField U_cell(cell_grid); U_cell = Zero();
//     localCopyRegion(U, U_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
//     LatticeLorentzScalar cell_mask = get_cell_mask(cell_grid);
//     U_cell = U_cell * cell_mask;
//
//
//     int n = 6; // FIXME: inner cell must be 6.6.6.6
// // # `lost` is the number of plaqs that are 0 because links perpendicular to the boundaries are set to 0
//     int lost = 4 * (n-1)*(n-1)*(n-1) * 3 + 6 * (n-1)*(n-1) * 5 + 4 * (n-1) * 6 + 6;
//
//     int total = n*n*n*n * 6;  // total number of size;
//     double multiplier = total / double(total - lost);  // Average plaq should be multiplied by this number 
//     std::cout << "Average Plaquette of U_cell: "<< WilsonLoops<PeriodicGimplR>::avgPlaquette(U_cell) * multiplier << std::endl;
//   }


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



}
