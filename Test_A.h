#include <Grid/Grid.h>
#include <Grid/Eigen/unsupported/MatrixFunctions>

#include "Test_gauge_force.h"

namespace Grid {
namespace QCD {

inline iMatrix<Complex, 3> Log(iMatrix<Complex, 3> &U, iVector<Complex, 3> &last_log_evals, const bool last_log_evals_initialized) { // parameter cannot be const in order to initialize in_eigen
  using namespace Eigen;
  // p.s. RowMajor only affect initialization
  using Mat = Eigen::Matrix<std::complex<double>, 3, 3, RowMajor>; // !!! by default MatrixXcd is columns Major !!!
  using Vec = Eigen::Vector3cd;

  Map<Mat> U_e(&U(0,0));

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

// inline iMatrix<Complex, 3> Log(iMatrix<Complex, 3> &in) { // parameter cannot be const in order to initialize in_eigen
//   using namespace Eigen;
//   iMatrix<Complex, 3> rst;
//   Map<MatrixXcd> in_eigen(&in(0,0), 3, 3);
//   MatrixXcd rst_eigen = in_eigen.log();
//
//   std::copy(rst_eigen.data(), rst_eigen.data() + 9, &rst(0,0) );
//   return rst;
// }

// // Does not work for beta = 10; U is not close to 1 enough
// template<class vtype, int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
// strong_inline iMatrix<vtype, N> Log(const iMatrix<vtype, N> &arg, Integer Nlog = 12 )
// // strong_inline iMatrix<vtype,N> Log(const iMatrix<vtype,N> &arg, Integer Nlog = 24 )
// {
//   typedef iMatrix<vtype,N> mat;
//   mat unit(1.0);
//   mat argMinusOne = arg - unit;
//
//   mat argMinusOnePowerI = unit;
//   mat temp = zero;
//   for(int i=1; i<=Nlog; ++i) {
//     argMinusOnePowerI *= argMinusOne;
//     if(i%2==1) temp += argMinusOnePowerI * (1. / double(i));
//     else temp -= argMinusOnePowerI * (1. / double(i));
//   }
//   return temp;
//
// }

using LatticeUeval = Lattice<iVector<iScalar<iVector<vComplex, 3> >, 4>>;
using LatticeUevalSite = iVector<iScalar<iVector<Complex, 3> >, 4>;

LatticeGaugeField Log(const LatticeGaugeField &lat, LatticeUeval &last_log_evals, bool &last_log_evals_initialized) {
  LatticeGaugeField rst(lat._grid);

  parallel_for(int ss=0; ss<lat._grid->lSites(); ss++) {
    std::vector<int> lcoor;
    lat._grid->LocalIndexToLocalCoor(ss, lcoor);

    LatticeGaugeFieldSite m;  LatticeUevalSite log_eval;
    peekLocalSite(m, lat, lcoor); peekLocalSite(log_eval, last_log_evals, lcoor);
    for(int mu=0; mu<4; ++mu) {
      m(mu)() = Log( m(mu)(), log_eval(mu)(), last_log_evals_initialized);
    }
    pokeLocalSite(log_eval, last_log_evals, lcoor);
    pokeLocalSite(m, rst, lcoor);


    // if(lcoor[0]==0 && lcoor[1]==0 && lcoor[2]==0 && lcoor[3]==0) std::cout << "[0,0,0,0] log eigenvalues: " << log_eval << std::endl;
    // if(lcoor[0]==1 && lcoor[1]==0 && lcoor[2]==0 && lcoor[3]==0) std::cout << "[1,0,0,0] log eigenvalues: " << log_eval << std::endl;
  }

  last_log_evals_initialized = true;

  return rst;
}


//
// LatticeGaugeField Log(const LatticeGaugeField &lat, ) {
//   LatticeGaugeField rst(lat._grid);
//
//   parallel_for(int ss=0; ss<lat._grid->lSites(); ss++) {
//     std::vector<int> lcoor;
//     lat._grid->LocalIndexToLocalCoor(ss, lcoor);
//
//     typename LatticeGaugeField::vector_object::scalar_object m;
//     peekLocalSite(m, lat, lcoor);
//     for(int mu=0; mu<4; ++mu) {
//       m(mu)() = Log( m(mu)() );
//     }
//     pokeLocalSite(m, rst, lcoor);
//   }
//
//   return rst;
// }

// LatticeColourMatrix Log(const LatticeColourMatrix &lat) {
//   LatticeColourMatrix rst(lat._grid);
//
//   parallel_for(int ss=0; ss<lat._grid->lSites(); ss++) {
//     std::vector<int> lcoor;
//     lat._grid->LocalIndexToLocalCoor(ss, lcoor);
//
//     typename LatticeColourMatrix::vector_object::scalar_object m;
//     peekLocalSite(m, lat, lcoor);
//     m()() = Log( m()() );
//     pokeLocalSite(m, rst, lcoor);
//   }
//   // parallel_for(int ss=0; ss<lat._grid->oSites(); ss++) {
//   //     rst[ss]()() = Log(lat[ss]()());
//   // }
//
//   return rst;
// }
//



void measure_A(const LatticeGaugeField &U, const std::vector<std::vector<int>> &coors, LatticeUeval &last_log_evals, bool &last_log_evals_initialized) {


    // std::cout << "before logrithm" << std::endl;
    // print_grid_field_site(U, {{0,0,0,0}});
    // print_grid_field_site(U, {{1,0,0,0}});


  LatticeGaugeField An(U._grid);
  An = timesMinusI(Log(U, last_log_evals, last_log_evals_initialized));

    // std::cout << "after logrithm" << std::endl;
    // LatticeGaugeField tmp = timesI(An);
    // print_grid_field_site(tmp, {{0,0,0,0}});
    // print_grid_field_site(tmp, {{1,0,0,0}});



  FFT theFFT((Grid::GridCartesian *)An._grid);
  theFFT.FFT_all_dim(An, An, FFT::forward);

  double vol = 1.0;
  for(int d=0; d<4; ++d)   vol = vol * U._grid->_fdimensions[d];
  An = An * (1. / std::sqrt(vol));

  // // FIXME: for test
  // std::cout << "after Fourier transform" << std::endl;
  // print_grid_field_site(An, {0,0,0,0});
  // print_grid_field_site(An, {1,0,0,0});

  Lattice<iVector<iScalar<iVector<vComplex, 8 > >, 4> > alg(U._grid);

  for (int mu = 0; mu < Nd; mu++)
  {
    LatticeColourMatrix An_mu(An._grid);
    SU3::LatticeAlgebraVector alg_mu(U._grid);

    An_mu = peekLorentz(An, mu);
    SU3::projectOnAlgebra(alg_mu, An_mu); // returns -2 i Tr(A t^a) = -i A^a; In Grid Tr(t^a t^b) = 1/2 \delta_{a,b}
    alg_mu = timesI(alg_mu);
    pokeLorentz(alg, alg_mu, mu);
  }
  // std::cout << alg << std::endl;

  for(const std::vector<int> &coor : coors) {
    // std::cout << "coor=" << coor << ": ";
    std::cout << "coor=";
    print_grid_field_site(alg, coor);
  }
}



// void measure_A(const LatticeGaugeField &U, const std::vector<std::vector<int>> &coors, bool log=true) {
//
//   LatticeGaugeField An(U._grid);
//   if(log) {                    // If input is U;
//     An = timesMinusI(Log(U));
//     // std::cout << "before logrithm" << std::endl;
//     // print_grid_field_site(U, {{0,0,0,0}});
//     // print_grid_field_site(U, {{1,0,0,0}});
//     // std::cout << "after logrithm" << std::endl;
//     // LatticeGaugeField tmp = timesI(An);
//     // print_grid_field_site(tmp, {, {1,0,0,0}});
//   }
//   else An = timesMinusI(U);   // If input is P;  In Grid, Mom is iP
//
//   An = An - trace(An) * (1. / 3.0);    // Remove trace; logrithm is ambiguous
//
//   FFT theFFT((Grid::GridCartesian *)An._grid);
//   // std::cout << "before Fourier transform" << std::endl;
//   // print_grid_field_site(An, {0,0,0,0});
//   // print_grid_field_site(An, {1,0,0,0});
//   theFFT.FFT_all_dim(An, An, FFT::forward);
//
//   double vol = 1.0;
//   for(int d=0; d<4; ++d)   vol = vol * U._grid->_fdimensions[d];
//   An = An * (1. / std::sqrt(vol));
//
//   // // FIXME: for test
//   // std::cout << "after Fourier transform" << std::endl;
//   // print_grid_field_site(An, {0,0,0,0});
//   // print_grid_field_site(An, {1,0,0,0});
//
//   Lattice<iVector<iScalar<iVector<vComplex, 8 > >, 4> > alg(U._grid);
//
//   for (int mu = 0; mu < Nd; mu++)
//   {
//     LatticeColourMatrix An_mu(An._grid);
//     SU3::LatticeAlgebraVector alg_mu(U._grid);
//
//     An_mu = peekLorentz(An, mu);
//     SU3::projectOnAlgebra(alg_mu, An_mu); // returns -2 i Tr(A t^a) = -i A^a; In Grid Tr(t^a t^b) = 1/2 \delta_{a,b}
//     alg_mu = timesI(alg_mu);
//     pokeLorentz(alg, alg_mu, mu);
//   }
//   // std::cout << alg << std::endl;
//
//   for(const std::vector<int> &coor : coors) {
//     // std::cout << "coor=" << coor << ": ";
//     std::cout << "coor=";
//     print_grid_field_site(alg, coor);
//   }
// }



void set_zero_mode_to_zero(LatticeGaugeField &P) {

  double V = 1.0;
  for(int d=0; d<4; ++d)   V = V * P._grid->_fdimensions[d];
  // std::cout << "V: " << V << std::endl;
  
  LatticeColourMatrix P_mu(P._grid);
  for(int mu=0; mu<4; mu++) {
    P_mu = peekLorentz(P, mu);

    typename LatticeColourMatrix::vector_object::scalar_object tmp;
    tmp = sum(P_mu);
    tmp = tmp * (1. / V); 
    
    // P_mu(n) = P_mu(n) - (1./V) * \sum_n' P_\mu(n');
    parallel_for(int ss=0; ss<P._grid->lSites(); ss++) {
      std::vector<int> lcoor;
      P._grid->LocalIndexToLocalCoor(ss, lcoor);

      typename LatticeColourMatrix::vector_object::scalar_object m;
      peekLocalSite(m, P_mu, lcoor);
      m = m - tmp;
      pokeLocalSite(m, P_mu, lcoor);
    }

    pokeLorentz(P, P_mu, mu);
  }
}


// void set_zero_mode_U_to_zero(LatticeGaugeField &U) {
//   std::cout << "Explicitly setting A(k=0) to 0" << std::endl;
//   using LatticeColourMatrixSite = typename LatticeColourMatrix::vector_object::scalar_object;
//
//   double V = 1.0;
//   for(int d=0; d<4; ++d)   V = V * U._grid->_fdimensions[d];
//   
//   LatticeColourMatrix U_mu(U._grid);
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






}}
