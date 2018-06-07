namespace Grid {
namespace QCD {

void GF_generate_P(LatticeGaugeField& P, GridParallelRNG& pRNG, const Momenta_k &KK)
{
  // step 1. generate Gaussian distributed matrix P0.
  LatticeGaugeField P0(P._grid);
  LatticeColourMatrix P0mu(P._grid);
  SU3::LatticeAlgebraVector P0_a(P._grid);   //Lattice<iScalar<iScalar<iVector<vComplex, SU3::AdjointDimension> > > >
  for (int mu = 0; mu < Nd; mu++) {
    gaussian(pRNG, P0_a); //P0_a's real and imaginary part both ~ N(0,1)
    SU3::FundamentalLieAlgebraMatrix(P0_a, P0mu);
    pokeLorentz(P0, P0mu, mu);
  }

  // step 2. transform P0
  LatticeColourMatrix sinKExpDotP0(P._grid);
  sinKExpDotP0 = KK.sinKPsExpDotP_func(P0);

  LatticeGaugeField newP(P._grid);
  LatticeColourMatrix Pmu(P._grid);
  for(int mu=0;mu<4;mu++){
  	Pmu = KK.sinKNgExp[mu] * KK.Ck_SqrtInvD * sinKExpDotP0 + KK.sinKEpsilon * peekLorentz(P0, mu);
  	pokeLorentz(newP, Pmu, mu);
  }

  //step 3. Impose reality condition $P_\mu(k) = \frac{1}{\sqrt{2}} (P_\mu(k) + P^\dagger_\mu(-k))$
  //generate P(-k) //FIXME: inefficient
  LatticeGaugeField newP_Minus(P._grid);

  std::vector<int> gdims = P._grid->_gdimensions;
  LorentzColourMatrix m;
  std::vector<int> coor;
  parallel_for(int ss=0;ss<P._grid->gSites();ss++){
    P._grid->GlobalIndexToGlobalCoor(ss, coor);
    peekSite(m, newP, coor);
    coor = gdims - coor; //zyd: do not need to modulo L. This is done in pokeSite.
    pokeSite(m, newP_Minus, coor);
  }

  // P_\mu(k) = \frac{1}{\sqrt{2}} (P_\mu(k) + P^\dagger_\mu(-k))
  LatticeGaugeField Pk(P._grid);
  Pk = newP + adj(newP_Minus);
  // Pk = Pk * (1.0/sqrt(2));
  Pk = Pk * 0.5;  //???? why

  //step 4.
  FFT theFFT((Grid::GridCartesian *)P._grid);
  theFFT.FFT_all_dim(P, Pk, FFT::backward);
  P = P * std::sqrt(KK.vol); //????is this appropriate

  // if(!isHermitian(P)) std::cout << "GF_generate_P:-----------yes hermitian------------" << std::endl;
  // if(!isHermitian(P)) std::cout << "GF_generate_P:-----------Not hermitian------------" << std::endl;
  // else std::cout << "GF_generate_P:-----------Yes hermitian------------" << std::endl;

  //In function GaussianFundamentalLieAlgebraMatrix, generated P is anti-hermitian; i * gaussian hermitian matrices
  P = timesI(P);

}


}}