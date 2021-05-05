#include <Grid/Grid.h>
#include <Grid/Eigen/unsupported/MatrixFunctions>

namespace Grid {

inline iMatrix<Complex, 3> Log(iMatrix<Complex, 3> &U, iVector<Complex, 3> &last_log_evals, const bool last_log_evals_initialized) { // parameter cannot be const in order to initialize in_eigen
  using namespace Eigen;
  // p.s. RowMajor only affect initialization
  using Mat = Eigen::Matrix<std::complex<double>, 3, 3, RowMajor>; // !!! by default MatrixXcd is columns Major !!!
  using Vec = Eigen::Vector3cd;

  // Grid::complex -> c++ complex
#if defined(GRID_CUDA) || defined(GRID_HIP)
  std::complex<double> U_tmp[3][3];
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      U_tmp[i][j] = std::complex<double>(U(i, j));
    }
  }


  Map<Mat> U_e(&U_tmp[0][0]);
  // Map<Mat> U_e(&U(0,0));
#else
  Map<Mat> U_e(&U(0,0));

#endif



  ComplexEigenSolver<Mat> es(U_e);
  Vec evals = es.eigenvalues();
  Mat evecs = es.eigenvectors();

  Vec log_evals = evals.array().log(); // Must add array() to do element-wise log; otherwise will call matrix log function

  // compare eigenvalues with the eigenvalues last step.
  for(int i=0; i<3; ++i) {
    if(last_log_evals_initialized) {   // if not first step
      // while(log_evals[i].imag() - last_log_evals(i).imag() >= M_PI) log_evals[i].imag() -= 2 * M_PI;
      while(log_evals[i].imag() - last_log_evals(i).imag() >= M_PI) log_evals[i].imag(log_evals[i].imag() - 2 * M_PI) ;
      while(log_evals[i].imag() - last_log_evals(i).imag() <= -M_PI) log_evals[i].imag(log_evals[i].imag() + 2 * M_PI) ;
      // while(log_evals[i].imag() - last_log_evals(i).imag() <= M_PI) log_evals[i].imag() += 2 * M_PI;
    }
    last_log_evals(i) = log_evals[i];
  }

  Mat rst_eigen = evecs * log_evals.asDiagonal() * evecs.adjoint();

  iMatrix<Complex, 3> rst;
  std::copy(rst_eigen.data(), rst_eigen.data() + 9, &rst(0,0) );
  return rst;
}


using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
using LatticeUevalSite = iVector<iScalar<iVector<Complex, 3> >, 4>;

LatticeGaugeField Log(const LatticeGaugeField &lat, LatticeUeval &last_log_evals, bool &last_log_evals_initialized) {
  LatticeGaugeField rst(lat.Grid());

  autoView(last_log_evals_v_read, last_log_evals, CpuRead);
  autoView(lat_v, lat, CpuRead);
  autoView(last_log_evals_v, last_log_evals, CpuWrite);
  autoView(rst_v, rst, CpuWrite);

  // parallel_for(int ss=0; ss<lat.Grid()->lSites(); ss++) {
  thread_for(ss, lat.Grid()->lSites(), {

    Coordinate lcoor;
    lat.Grid()->LocalIndexToLocalCoor(ss, lcoor);

    LatticeGaugeFieldSite m;  LatticeUevalSite log_eval;
    peekLocalSite(m, lat_v, lcoor); peekLocalSite(log_eval, last_log_evals_v_read, lcoor);
    for(int mu=0; mu<4; ++mu) {
      m(mu)() = Log( m(mu)(), log_eval(mu)(), last_log_evals_initialized);
    }

    pokeLocalSite(log_eval, last_log_evals_v, lcoor);
    pokeLocalSite(m, rst_v, lcoor);

    // if(lcoor[0]==0 && lcoor[1]==0 && lcoor[2]==0 && lcoor[3]==0) std::cout << "[0,0,0,0] log eigenvalues: " << log_eval << std::endl;
    // if(lcoor[0]==1 && lcoor[1]==0 && lcoor[2]==0 && lcoor[3]==0) std::cout << "[1,0,0,0] log eigenvalues: " << log_eval << std::endl;
  });

  last_log_evals_initialized = true;

  return rst;
}



void measure_A(const LatticeGaugeField &U, const std::vector<std::vector<int>> &coors, LatticeUeval &last_log_evals, bool &last_log_evals_initialized) {


  LatticeGaugeField An(U.Grid());
  An = timesMinusI(Log(U, last_log_evals, last_log_evals_initialized));

  FFT theFFT((Grid::GridCartesian *)An.Grid());
  theFFT.FFT_all_dim(An, An, FFT::forward);

  double vol = 1.0;
  for(int d=0; d<4; ++d)   vol = vol * U.Grid()->_fdimensions[d];
  An = An * (1. / std::sqrt(vol));

  Lattice<iVector<iScalar<iVector<vComplex, 8 > >, 4> > alg(U.Grid());

  for (int mu = 0; mu < Nd; mu++)
  {
    LatticeColourMatrix An_mu(An.Grid());
    SU3::LatticeAlgebraVector alg_mu(U.Grid());

    An_mu = peekLorentz(An, mu);
    SU3::projectOnAlgebra(alg_mu, An_mu); // returns -2 i Tr(A t^a) = -i A^a; In Grid Tr(t^a t^b) = 1/2 \delta_{a,b}
    alg_mu = timesI(alg_mu);
    pokeLorentz(alg, alg_mu, mu);
  }

  for(const std::vector<int> &coor : coors) {
    // std::cout << "coor=" << coor << ": ";
    std::cout << "coor=";
    print_grid_field_site(alg, coor);
  }
  // std::cout << "end of measure_A function" << std::endl;
}



void set_zero_mode_to_zero(LatticeGaugeField &P) {

  double V = 1.0;
  for(int d=0; d<4; ++d)   V = V * P.Grid()->_fdimensions[d];
  // std::cout << "V: " << V << std::endl;
  
  LatticeColourMatrix P_mu(P.Grid());
  for(int mu=0; mu<4; mu++) {
    P_mu = peekLorentz(P, mu);

    typename LatticeColourMatrix::vector_object::scalar_object tmp;
    tmp = sum(P_mu);
    tmp = tmp * (1. / V); 
    
    // P_mu(n) = P_mu(n) - (1./V) * \sum_n' P_\mu(n');

    autoView(P_mu_v_write, P_mu, CpuWrite);
    thread_for(ss, P.Grid()->lSites(), {
      Coordinate lcoor;
      P.Grid()->LocalIndexToLocalCoor(ss, lcoor);

      typename LatticeColourMatrix::vector_object::scalar_object m;
      peekLocalSite(m, P_mu.View(CpuRead), lcoor);
      m = m - tmp;
      pokeLocalSite(m, P_mu_v_write, lcoor);
    });

    pokeLorentz(P, P_mu, mu);
  }
}


// void set_zero_mode_U_to_zero(LatticeGaugeField &U) {
//   std::cout << "Explicitly setting A(k=0) to 0" << std::endl;
//   using LatticeColourMatrixSite = typename LatticeColourMatrix::vector_object::scalar_object;
//
//   double V = 1.0;
//   for(int d=0; d<4; ++d)   V = V * U.Grid()->_fdimensions[d];
//   
//   LatticeColourMatrix U_mu(U.Grid());
//   for(int mu=0; mu<4; mu++) {
//     U_mu = peekLorentz(U, mu);
//
//     LatticeColourMatrix log_Umu = Log(U_mu);
//     LatticeColourMatrixSite tmp;
//     tmp = sum(log_Umu);
//     tmp = tmp * (-1. / V); 
//     tmp = Ta(tmp); // project onto su(3) lie algebra // remove the trace and make it antisymmetric
//
//     LatticeColourMatrixSite exp_mat = ProjectOnGroup(Exponentiate(tmp, 1));
//
//     U_mu = exp_mat * U_mu;
//     
//     pokeLorentz(U, U_mu, mu);
//     
//   }
// }






}
