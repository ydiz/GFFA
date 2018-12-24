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
  std::vector<int> coor0000;
  std::vector<LatticeComplex> sinK; // sin(k/2)
  std::vector<LatticeComplex> sinKNgExp; // sin(k/2) exp(-ik/2)
  std::vector<LatticeComplex> sinKPsExp; // sin(k/2) exp(ik/2)
  LatticeComplex sinKNormSquare;
  LatticeComplex sinKNorm;
  LatticeComplex sinKEpsilonSquare; //sin(k/2)^2+\epsilon^2
  LatticeComplex sinKEpsilon; //\sqrt{sin(k/2)^2+\epsilon^2}
  LatticeComplex Ck_D; //  (\frac{1}{M^2} - \frac{1}{sin(k/2)^2+\epsilon^2})\frac{1}{sin(k/2)^2}
  LatticeComplex Ck_SqrtInvD; //(M-\sqrt{sin(k/2)^2+\epsilon^2})\frac{1}{sin(k/2)^2}


  Momenta_k(GridBase* grid, Real _M, Real _epsi, bool _newHp): k(4, grid), sinK(4, grid), sinKNgExp(4, grid), sinKPsExp(4, grid), coor0000(4), sinKNormSquare(grid),
  sinKNorm(grid), sinKEpsilonSquare(grid), sinKEpsilon(grid),  Ck_D(grid), Ck_SqrtInvD(grid), M(_M), newHp(_newHp)
  {
    // generate coor0000, k, and sinK
    LatticeComplex    coor(grid);
    //LatticeComplex    coor2(grid);
    LatticeComplex    tt(grid);
    for(int mu=0;mu<4;mu++){
      coor0000[mu] =  grid->_fdimensions[mu] / 2;
      RealD TwoPiL =  M_PI * 2.0/ grid->_fdimensions[mu];
      LatticeCoordinate(coor,mu);
	  assign_lc(tt,  std::complex<double>(coor0000[mu],0));
	  //coor2 = coor - tt;//std::complex<double>(coor0000[mu],0);
      k[mu] = TwoPiL * (coor - tt);//(coor - std::complex<double>(coor0000[mu],0));
      sinK[mu] = sin(k[mu] * 0.5);
      sinKPsExp[mu] = exp(timesI(k[mu]*0.5)) * sinK[mu];
      sinKNgExp[mu] = adj(sinKPsExp[mu]);
    }

    //generate epsi;
    epsi = _epsi;

    iScalar<iScalar<iScalar<vComplex> >> epsilon; //regulator; avoid 0 in denominator k.
    epsilon = epsi;

    //generate sinKNormSquare
    sinKNormSquare = sinK[0]*sinK[0] +sinK[1]*sinK[1] + sinK[2]*sinK[2] +sinK[3]*sinK[3];
    sinKEpsilonSquare = sinKNormSquare + epsilon*epsilon;
    // sinKEpsilonSquare = epsilon*epsilon;

    //generate sinKEpsilon
    sinKEpsilon = sqrt(sinKEpsilonSquare);

    //generate one
    one = 1.0;

    //generate Ck_D.
    TComplex zero;
    zero = 0.0;
    Ck_D = one / (sinKNormSquare * (M * M)) - one / (sinKEpsilonSquare * sinKNormSquare);
    pokeSite(zero, Ck_D, coor0000);

    //generate Ck_SqrtInvD
    Ck_SqrtInvD = (M-sinKEpsilon)/sinKNormSquare;
    pokeSite(zero, Ck_SqrtInvD, coor0000);

    //volume vol
    vol = 1.0;
    for(int d=0; d<Nd; ++d)   vol = vol * grid->_fdimensions[d];

    // std::cout << __func__ << ": M = " << M << std::endl;
    // std::cout << __func__ << ": epsi = " << epsi << std::endl;
    // std::cout << __func__ << ": newHp = " << newHp << std::endl;

  }

  LatticeColourMatrix sinKPsExpDotP_func(const LatticeGaugeField &P) const
  {
    LatticeColourMatrix kDotP(P._grid);
    kDotP = zero;
    for(int mu=0;mu<4;mu++){
      kDotP += sinKPsExp[mu] * PeekIndex<LorentzIndex>(P, mu);
    }
    return kDotP;
  }


};



}
}
