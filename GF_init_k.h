namespace Grid{
namespace QCD{

class Momenta_k{

public:
  Real M;
  Real vol;
  vTComplex one;
  Real epsi;
  bool newHp;
  std::vector<LatticeComplex> k;//(4, grid); //cannot initialzie here
  std::vector<LatticeComplex> sinK; // sin(k/2)
  std::vector<LatticeComplex> sinKNgExp; // sin(k/2) exp(-ik/2)
  std::vector<LatticeComplex> sinKPsExp; // sin(k/2) exp(ik/2)
  LatticeComplex sinKNormSquare;
  LatticeComplex sinKNorm;
  LatticeComplex FourSinKSquareEpsilon; //sin(k/2)^2+\epsilon^2
  LatticeComplex SqrtFourSinKSquareEpsilon; //sin(k/2)^2+\epsilon^2
  LatticeComplex Ck_D; //  (\frac{1}{M^2} - \frac{1}{sin(k/2)^2+\epsilon^2})\frac{1}{sin(k/2)^2}
  LatticeComplex Ck_SqrtInvD; //(M-\sqrt{sin(k/2)^2+\epsilon^2})\frac{1}{sin(k/2)^2}


  Momenta_k(GridBase* grid, Real _M, Real _epsi, bool _newHp): k(4, grid), sinK(4, grid), sinKNgExp(4, grid), sinKPsExp(4, grid), sinKNormSquare(grid), sinKNorm(grid), FourSinKSquareEpsilon(grid), SqrtFourSinKSquareEpsilon(grid), Ck_D(grid), Ck_SqrtInvD(grid), M(_M), newHp(_newHp)
  {
    // generate k, and sinK
    LatticeComplex    coor(grid);
    LatticeComplex    tt(grid);
    for(int mu=0;mu<4;mu++){
      RealD TwoPiL =  M_PI * 2.0/ grid->_fdimensions[mu];
      LatticeCoordinate(coor,mu);
      k[mu] = TwoPiL * coor;
      sinK[mu] = sin(k[mu] * 0.5);
      sinKPsExp[mu] = exp(timesI(k[mu]*0.5)) * sinK[mu];
      sinKNgExp[mu] = adj(sinKPsExp[mu]);
    }

    //generate epsi;
    epsi = _epsi;


    //generate sinKNormSquare
    sinKNormSquare = sinK[0]*sinK[0] +sinK[1]*sinK[1] + sinK[2]*sinK[2] +sinK[3]*sinK[3];
    FourSinKSquareEpsilon = 4. * sinKNormSquare; // + epsilon*epsilon;

    iScalar<iScalar<iScalar<Complex> >> epsilon; //regulator; avoid 0 in denominator k.
    epsilon = epsi;
    pokeSite(epsilon * epsilon, FourSinKSquareEpsilon, Coordinate({0,0,0,0}));  // Give zero mode mass epsilon

    

    //generate sinKEpsilon
    SqrtFourSinKSquareEpsilon = sqrt(FourSinKSquareEpsilon);

    //generate one
    one = 1.0;

    //generate Ck_D.
    TComplex zero;
    zero = 0.0;
    Ck_D = one / (sinKNormSquare * (M * M)) - one / (FourSinKSquareEpsilon * sinKNormSquare);
    pokeSite(zero, Ck_D, Coordinate({0,0,0,0}));

    //generate Ck_SqrtInvD
    Ck_SqrtInvD = (M-SqrtFourSinKSquareEpsilon)/sinKNormSquare;
    // pokeSite(zero, Ck_SqrtInvD, coor0000);
    pokeSite(zero, Ck_SqrtInvD, Coordinate({0,0,0,0}));

    //volume vol
    vol = 1.0;
    for(int d=0; d<Nd; ++d)   vol = vol * grid->_fdimensions[d];

    // std::cout << __func__ << ": M = " << M << std::endl;
    // std::cout << __func__ << ": epsi = " << epsi << std::endl;
    // std::cout << __func__ << ": newHp = " << newHp << std::endl;

  }

  LatticeColourMatrix sinKNgExpDotP_func(const LatticeGaugeField &P) const
  {
    LatticeColourMatrix kDotP(P.Grid());
    kDotP = Zero();
    for(int mu=0;mu<4;mu++){
      kDotP += sinKNgExp[mu] * PeekIndex<LorentzIndex>(P, mu);
    }
    return kDotP;
  }


};



}
}
