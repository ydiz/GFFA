namespace Grid{
  namespace QCD{
    //take from Test_GaugeAction.cc
bool metropolis_test(const RealD Delta, GridSerialRNG &sRNG) {
  RealD rand;
  if(Delta <=0.0) return true;

  random(sRNG,rand);
  if(rand <= std::exp(-Delta))
    return true;
  else
    return false;
}

  // modified from functions in SUn.h
template <typename MatrixType>
static void taExp(const MatrixType &x, MatrixType &ex) {
  typedef typename MatrixType::scalar_type ComplexType;

  MatrixType xn;
  RealD nfac = 1.0;

  xn = x;
  ex = xn + ComplexType(1.0);  // 1+x

  // Do a 12th order exponentiation
  for (int i = 2; i <= 12; ++i) {
    nfac = nfac / RealD(i);  // 1/2, 1/2.3 ...
    xn = xn * x;             // x2, x3,x4....
    ex = ex + xn * nfac;     // x2/2!, x3/3!....
  }
}

// template <typename MatrixType>
static void LieRandomize_serial(GridSerialRNG &sRNG, ColourMatrix &out,
                         double scale = 1.0) {
  Complex ca;
  ColourMatrix lie;
  ColourMatrix la;
  Complex ci(0.0, scale);
  Complex cone(1.0, 0.0);
  ColourMatrix  ta;

  lie = zero;
  for (int a = 0; a < 8; a++) {
    random(sRNG, ca);
    ca = (ca + conjugate(ca)) * 0.5;
    ca = ca - 0.5;
    SU<3>::generator(a, ta);

    // ci * ca * ta;
    la = ci * ca * ta;

    lie = lie + la;  // e^{i la ta}
  }
  taExp(lie, out);
}

}
}
