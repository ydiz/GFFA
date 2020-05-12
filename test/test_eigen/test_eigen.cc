#include <Grid/Grid.h>
// #include <Grid/Eigen/unsupported/src/MatrixFunctions/MatrixLogarithm.h>
#include <Grid/Eigen/unsupported/MatrixFunctions>



using namespace std;
using namespace Eigen;

int main() {

  // p.s. RowMajor only affect initialization
  using Mat = Matrix<std::complex<double>, 3, 3, RowMajor>; // !!! by default MatrixXcd is columns Major !!!
  using Vec = Vector3cd;

  const std::complex<double> j(0, 1);
  complex<double> a[3][3] = {{-0.41707526 + 0.18242808*j, -0.55927303 + -0.32881545*j, 0.16997965 + -0.58563573*j}, {-0.7906834 + -0.00092663*j, -0.08918532 + 0.41383915*j, 0.07165142 + 0.43642653*j}, {-0.12377057 + 0.39021528*j, -0.02543847 + -0.63184658*j, -0.4802048 + 0.44937625*j}};


  Map<Mat> U(&a[0][0]); // If not set RowMajor, U will be the transpose of a.

  std::cout << U << std::endl;

  ComplexEigenSolver<Mat> es(U);

  Vec evals = es.eigenvalues();
  Mat evecs = es.eigenvectors();
  std::cout << evals << std::endl;

  Vec evals_log = evals.array().log();
  std::cout << evals_log << std::endl;

  std::cout << evecs * evals_log.asDiagonal() * evecs.adjoint() << std::endl;

};

