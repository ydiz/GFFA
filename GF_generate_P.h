namespace Grid {
namespace QCD {

void GF_generate_P(LatticeGaugeField& P, GridParallelRNG& pRNG, const Momenta_k &KK)
{
  GridBase *full_grid = P.Grid();
  GridBase *cell_grid = KK.Grid();
  Coordinate cell_size  = cell_grid->_fdimensions;

  // step 1. generate Gaussian distributed matrix P0.
  LatticeGaugeField P0(P.Grid());
  LatticeColourMatrix P0mu(P.Grid());
  SU3::LatticeAlgebraVector P0_a(P.Grid());   //Lattice<iScalar<iScalar<iVector<vComplex, SU3::AdjointDimension> > > >
  for(int mu = 0; mu < Nd; mu++) {
    gaussian(pRNG, P0_a); //P0_a's real and imaginary part both ~ N(0,1)
    SU3::FundamentalLieAlgebraMatrix(P0_a, P0mu);
    pokeLorentz(P0, P0mu, mu);
  }

  LatticeGaugeField P0_cell(cell_grid);
  localCopyRegion(P0, P0_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

  // step 2. transform P0
  LatticeColourMatrix sinKExpDotP0 = KK.sinKNgExpDotP_func(P0_cell);

  LatticeGaugeField newP(cell_grid);
  LatticeColourMatrix Pmu(cell_grid);
  for(int mu=0;mu<4;mu++){
  	Pmu = KK.sinKPsExp[mu] * KK.Ck_SqrtInvD * sinKExpDotP0 + KK.SqrtFourSinKSquareEpsilon * peekLorentz(P0_cell, mu);
  	pokeLorentz(newP, Pmu, mu);
  }

  // // P_\mu(k) = \frac{1}{\sqrt{2}} (P_\mu(k) + P^\dagger_\mu(-k))
  LatticeGaugeField Pk(cell_grid);
  // Pk = newP + adj(newP_Minus);
  Pk = newP;  // Do not need to impose reality condition on Pk; it is equivalent to do it directly on P(x) // sum_k P^\dagger_\mu(-k) e^{ikx} = sum_k P^\dagger_\mu(k) e^{-ikx} = ( sum_k P_\mu(k) e^{ikx}  )^\dagger

  //step 4.
  LatticeGaugeField P_cell(cell_grid);
  FFT theFFT((Grid::GridCartesian *)cell_grid);
  theFFT.FFT_all_dim(P_cell, Pk, FFT::backward);
  P_cell = P_cell * std::sqrt(KK.vol); 

  P_cell = (P_cell + adj(P_cell)) * 0.5;  // P must be hermitian // Must multiply it by 0.5 
  P_cell = P_cell * sqrt(HMC_MOMENTUM_DENOMINATOR);   // IMPORTANT: after Grid changed its convention since Dec 19 2018, must multiply the initial momentum by sqrt(2).   See generate_momenta function in Grid/qcd/action/gauge/GaugeImplTypes.h
  // if(!isHermitian(P)) std::cout << "GF_generate_P:-----------Not hermitian------------" << std::endl;

  //In GaussianFundamentalLieAlgebraMatrix, P is anti-hermitian; i * gaussian hermitian matrices
  P_cell = timesI(P_cell);


  //  insert P_cell into P

  P = Zero();
  localCopyRegion(P_cell, P, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);
}


}}
