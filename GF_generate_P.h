namespace Grid {
namespace QCD {

void GF_generate_P(LatticeGaugeField& P, GridParallelRNG& pRNG, const Momenta_k &KK)
{
  // step 1. generate Gaussian distributed matrix P0.
  LatticeGaugeField P0(P._grid);
  LatticeColourMatrix P0mu(P._grid);
  SU3::LatticeAlgebraVector P0_a(P._grid);   //Lattice<iScalar<iScalar<iVector<vComplex, SU3::AdjointDimension> > > >
  for(int mu = 0; mu < Nd; mu++) {
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

// If use parallel_for, it will generate an error; do not know why
  for(int node=0; node<P._grid->_Nprocessors; ++node){
	for(int ss=0; ss<P._grid->lSites(); ss++){
		LorentzColourMatrix m;
		std::vector<int> lcoor(4);
		P._grid->LocalIndexToLocalCoor(ss, lcoor);
		std::vector<int> gcoor(4);
		std::vector<int> processor_coor;
		P._grid->ProcessorCoorFromRank(node, processor_coor); // to use peekSite, all nodes must put in the same coor
		P._grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
		//std::cout << processor_coor << lcoor << gcoor << std::endl;
		std::vector<int> new_gcoor(4);
		new_gcoor = gdims - gcoor;
		peekSite(m, newP, new_gcoor);
		//gcoor = gdims - gcoor; //zyd: do not need to modulo L. This is done in pokeSite.
		// std::cout << gcoor << std::endl;
		if(P._grid->ThisRank()==node) {
			pokeLocalSite(m, newP_Minus, lcoor);
		}
	}
	MPI_Barrier(P._grid->communicator_world);

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
