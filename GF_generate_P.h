namespace Grid {
namespace QCD {

void GF_generate_P(LatticeGaugeField& P, GridParallelRNG& pRNG, const Momenta_k &KK)
{
  
  // step 1. generate Gaussian distributed matrix P0.
  LatticeGaugeField P0(P.Grid());
  LatticeColourMatrix P0mu(P.Grid());
  SU3::LatticeAlgebraVector P0_a(P.Grid());   //Lattice<iScalar<iScalar<iVector<vComplex, SU3::AdjointDimension> > > >
  for(int mu = 0; mu < Nd; mu++) {
    gaussian(pRNG, P0_a); //P0_a's real and imaginary part both ~ N(0,1)
    SU3::FundamentalLieAlgebraMatrix(P0_a, P0mu);
    pokeLorentz(P0, P0mu, mu);
  }

  // step 2. transform P0
  LatticeColourMatrix sinKExpDotP0(P.Grid());
  sinKExpDotP0 = KK.sinKNgExpDotP_func(P0);

  LatticeGaugeField newP(P.Grid());
  LatticeColourMatrix Pmu(P.Grid());
  for(int mu=0;mu<4;mu++){
  	Pmu = KK.sinKPsExp[mu] * KK.Ck_SqrtInvD * sinKExpDotP0 + KK.SqrtFourSinKSquareEpsilon * peekLorentz(P0, mu);
  	pokeLorentz(newP, Pmu, mu);
  }

//   //step 3. Impose reality condition $P_\mu(k) = \frac{1}{\sqrt{2}} (P_\mu(k) + P^\dagger_\mu(-k))$
//   //generate P(-k) //FIXME: inefficient
//   LatticeGaugeField newP_Minus(P.Grid());
//
//   Coordinate gdims = P.Grid()->_gdimensions;
//
//   autoView(newP_Minus_v, newP_Minus, CpuWrite);
//
// // Cannot use parallel_for 
//   for(int node=0; node<P.Grid()->_Nprocessors; ++node){
// 	for(int ss=0; ss<P.Grid()->lSites(); ss++){
// 		LorentzColourMatrix m;
// 		// std::vector<int> lcoor(4);
//     Coordinate lcoor, gcoor, processor_coor, new_gcoor;
// 		P.Grid()->LocalIndexToLocalCoor(ss, lcoor);
// 		P.Grid()->ProcessorCoorFromRank(node, processor_coor); // to use peekSite, all nodes must put in the same coor
// 		P.Grid()->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
//
// 		new_gcoor = gdims - gcoor;
// 		peekSite(m, newP, new_gcoor);
//
// 		//gcoor = gdims - gcoor; // zyd: do not need to modulo L. This is done in pokeSite.
// 		if(P.Grid()->ThisRank()==node) {
// 			// pokeLocalSite(m, newP_Minus, lcoor);
// 			pokeLocalSite(m, newP_Minus_v, lcoor);
// 		}
// 	}
// #ifndef GRID_COMMS_NONE
// 	MPI_Barrier(P.Grid()->communicator_world);
// #endif
//   }
  //
  // // P_\mu(k) = \frac{1}{\sqrt{2}} (P_\mu(k) + P^\dagger_\mu(-k))
  LatticeGaugeField Pk(P.Grid());
  // Pk = newP + adj(newP_Minus);
  Pk = newP;  // Do not need to impose reality condition on Pk; it is equivalent to do it directly on P(x) // sum_k P^\dagger_\mu(-k) e^{ikx} = sum_k P^\dagger_\mu(k) e^{-ikx} = ( sum_k P_\mu(k) e^{ikx}  )^\dagger

  //step 4.
  FFT theFFT((Grid::GridCartesian *)P.Grid());
  theFFT.FFT_all_dim(P, Pk, FFT::backward);
  P = P * std::sqrt(KK.vol); 

  P = (P + adj(P)) * 0.5;  // P must be hermitian // Must multiply it by 0.5 
  P = P * sqrt(HMC_MOMENTUM_DENOMINATOR);   // IMPORTANT: after Grid changed its convention since Dec 19 2018, must multiply the initial momentum by sqrt(2).   See generate_momenta function in Grid/qcd/action/gauge/GaugeImplTypes.h
  // if(!isHermitian(P)) std::cout << "GF_generate_P:-----------Not hermitian------------" << std::endl;

  //In GaussianFundamentalLieAlgebraMatrix, P is anti-hermitian; i * gaussian hermitian matrices
  P = timesI(P);
}


}}
