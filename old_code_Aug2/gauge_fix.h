namespace Grid {
namespace QCD {

// FIXME: haven't changed g in FourierAccelSteepestDescentStep!!!!!!!!!!!!!
// added paramter g to SteepestDescentGaugeFix, SteepestDescentStep, and FourierAccelSteepestDescentStep
// added: update g as well after update U
template <class Gimpl>
class My_FourierAcceleratedGaugeFixer  : public Gimpl {
 public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Gimpl::GaugeLinkField GaugeMat;
  typedef typename Gimpl::GaugeField GaugeLorentz;

  static void GaugeLinkToLieAlgebraField(const std::vector<GaugeMat> &U,std::vector<GaugeMat> &A) {
    for(int mu=0;mu<Nd;mu++){
      Complex cmi(0.0,-1.0);
      A[mu] = Ta(U[mu]) * cmi;
    }
  }
  static void DmuAmu(const std::vector<GaugeMat> &A,GaugeMat &dmuAmu) {
    dmuAmu=zero;
    for(int mu=0;mu<Nd;mu++){
      dmuAmu = dmuAmu + A[mu] - Cshift(A[mu],mu,-1);
    }
  }
  static void SteepestDescentGaugeFix(GaugeLorentz &Umu,GaugeMat &g, Real & alpha,int maxiter,Real Omega_tol, Real Phi_tol,bool Fourier=false) {
    GridBase *grid = Umu._grid;

    Real org_plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
    Real org_link_trace=WilsonLoops<Gimpl>::linkTrace(Umu);
    Real old_trace = org_link_trace;
    Real trG;

    std::vector<GaugeMat> U(Nd,grid);
    GaugeMat dmuAmu(grid);

    for(int i=0;i<maxiter;i++){
      for(int mu=0;mu<Nd;mu++) U[mu]= PeekIndex<LorentzIndex>(Umu,mu);
      if ( Fourier==false ) {
	trG = SteepestDescentStep(U, g, alpha,dmuAmu);
      } else {
	trG = FourierAccelSteepestDescentStep(U, g, alpha,dmuAmu);
      }
      for(int mu=0;mu<Nd;mu++) PokeIndex<LorentzIndex>(Umu,U[mu],mu);
      // Monitor progress and convergence test
      // infrequently to minimise cost overhead
      if ( i %20 == 0 ) {
      	Real plaq      =WilsonLoops<Gimpl>::avgPlaquette(Umu);
      	Real link_trace=WilsonLoops<Gimpl>::linkTrace(Umu);

      	// if (Fourier)
      	//   std::cout << GridLogMessage << "Fourier Iteration "<<i<< " plaq= "<<plaq<< " dmuAmu " << norm2(dmuAmu)<< std::endl;
      	// else
      	//   std::cout << GridLogMessage << " Iteration "<<i<< " plaq= "<<plaq<< " dmuAmu " << norm2(dmuAmu)<< std::endl;

      	Real Phi  = 1.0 - old_trace / link_trace ;
      	Real Omega= 1.0 - trG;


      	// std::cout << GridLogMessage << " Iteration "<<i<< " Phi= "<<Phi<< " Omega= " << Omega<< " trG " << trG <<std::endl;
      	if ( (Omega < Omega_tol) && ( ::fabs(Phi) < Phi_tol) ) {
      	// if (  ::fabs(Phi) < Phi_tol ) {
      	  std::cout << GridLogMessage << "Converged ! "<<std::endl;
      	  return;
      	}

      	old_trace = link_trace;

      }
    }
  };
  static Real SteepestDescentStep(std::vector<GaugeMat> &U, GaugeMat &g, Real & alpha, GaugeMat & dmuAmu) {
    GridBase *grid = U[0]._grid;

    std::vector<GaugeMat> A(Nd,grid);
    GaugeMat g_tmp(grid);

    GaugeLinkToLieAlgebraField(U,A);
    ExpiAlphaDmuAmu(A,g_tmp,alpha,dmuAmu);


    Real vol = grid->gSites();
    Real trG = TensorRemove(sum(trace(g_tmp))).real()/vol/Nc;

    SU<Nc>::GaugeTransform(U,g_tmp);
    g = g_tmp * g;

    return trG;
  }

  static Real FourierAccelSteepestDescentStep(std::vector<GaugeMat> &U, GaugeMat &g, Real & alpha, GaugeMat & dmuAmu) {

    GridBase *grid = U[0]._grid;

    Real vol = grid->gSites();

    FFT theFFT((GridCartesian *)grid);

    LatticeComplex  Fp(grid);
    LatticeComplex  psq(grid); psq=zero;
    LatticeComplex  pmu(grid);
    LatticeComplex   one(grid); one = Complex(1.0,0.0);

    // GaugeMat g(grid);
    GaugeMat dmuAmu_p(grid);
    std::vector<GaugeMat> A(Nd,grid);

    GaugeLinkToLieAlgebraField(U,A);

    DmuAmu(A,dmuAmu);

    theFFT.FFT_all_dim(dmuAmu_p,dmuAmu,FFT::forward);

    //////////////////////////////////
    // Work out Fp = psq_max/ psq...
    //////////////////////////////////
    std::vector<int> latt_size = grid->GlobalDimensions();
    std::vector<int> coor(grid->_ndimension,0);
    for(int mu=0;mu<Nd;mu++) {

      Real TwoPiL =  M_PI * 2.0/ latt_size[mu];
      LatticeCoordinate(pmu,mu);
      pmu = TwoPiL * pmu ;
      psq = psq + 4.0*sin(pmu*0.5)*sin(pmu*0.5);
    }

    Complex psqMax(16.0);
    Fp =  psqMax*one/psq;

    /*
    static int once;
    if ( once == 0 ) {
      std::cout << " Fp " << Fp <<std::endl;
      once ++;
      }*/

    pokeSite(TComplex(1.0),Fp,coor);

    dmuAmu_p  = dmuAmu_p * Fp;

    theFFT.FFT_all_dim(dmuAmu,dmuAmu_p,FFT::backward);

    GaugeMat ciadmam(grid);
    Complex cialpha(0.0,-alpha);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);

    Real trG = TensorRemove(sum(trace(g))).real()/vol/Nc;

    SU<Nc>::GaugeTransform(U,g);

    return trG;
  }

  static void ExpiAlphaDmuAmu(const std::vector<GaugeMat> &A,GaugeMat &g,Real & alpha, GaugeMat &dmuAmu) {
    GridBase *grid = g._grid;
    Complex cialpha(0.0,-alpha);
    GaugeMat ciadmam(grid);
    DmuAmu(A,dmuAmu);
    ciadmam = dmuAmu*cialpha;
    SU<Nc>::taExp(ciadmam,g);
  }
};

}
}
