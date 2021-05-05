#include <cmath>

namespace Grid{

//calculate coordinate after shifting considering periodic boundary condition
template<class vobj> int shiftUp(const Lattice<vobj> &U, int x, int mu)
{
	std::vector<int> coor;
  Lexicographic::CoorFromIndex(coor, x, U.Grid()->_gdimensions);
	int gd = U.Grid()->_gdimensions[mu];
	if(coor[mu]==gd-1) return x - (gd-1)*U.Grid()->_ostride[mu];
	else return x+U.Grid()->_ostride[mu];
}

template<class vobj> int shiftDown(const Lattice<vobj> &U, int x, int mu)
{
  std::vector<int> coor;
  Lexicographic::CoorFromIndex(coor, x, U.Grid()->_gdimensions);
  int gd = U.Grid()->_gdimensions[mu];
  if(coor[mu]==0) return x + (gd-1)*U.Grid()->_ostride[mu];
  else return x-U.Grid()->_ostride[mu];
}

//didn't consider mpi.
//peek vColourMatrix
vColourMatrix peekvCMatrix(const LatticeGaugeField& U, int x, int mu)
{
  int nsimd = U.Grid()->Nsimd();
  vColourMatrix vM;

  int osites = U.Grid()->oSites();
  int gsites = U.Grid()->gSites();
  std::vector<int> gdims = U.Grid()->_gdimensions;

  std::vector<std::vector<int>> coor(nsimd);
  for(int idx=0; idx<nsimd; ++idx)
  {
    coor[idx].resize(Nd); //??redundancy??
    Lexicographic::CoorFromIndex(coor[idx], (x + idx*osites)%gsites, gdims);
  }

  std::vector<ColourMatrix> Mvec(nsimd);
  std::vector<LorentzColourMatrix> LMvec(nsimd);
  for(int idx=0; idx<nsimd; ++idx)
  {
    peekSite(LMvec[idx], U, coor[idx]);
    Mvec[idx] = peekIndex<LorentzIndex>(LMvec[idx], mu);
  }

  merge(vM, Mvec);
  return vM;
}

vColourMatrix peekvg(const LatticeColourMatrix& g, int x)
{
  int nsimd = g.Grid()->Nsimd();
  vColourMatrix vM;

  int osites = g.Grid()->oSites();
  int gsites = g.Grid()->gSites();
  std::vector<int> gdims = g.Grid()->_gdimensions;

  std::vector<std::vector<int>> coor(nsimd);
  for(int idx=0; idx<nsimd; ++idx)
  {
    coor[idx].resize(Nd);
    Lexicographic::CoorFromIndex(coor[idx], (x + idx*osites)%gsites, gdims);
  }

  std::vector<ColourMatrix> Mvec(nsimd);
  for(int idx=0; idx<nsimd; ++idx)
  {
    peekSite(Mvec[idx], g, coor[idx]);
  }

  merge(vM, Mvec);
  return vM;
}


void update_gx(int x, LatticeColourMatrix &g, const LatticeGaugeField &U, GridSerialRNG &sRNG, RealD betaMM, RealD scale, int multi_hit)
{
  int nsimd = g.Grid()->Nsimd();

  vColourMatrix gx_new;
  vColourMatrix vM; //random matrix used to update Uxmu
  std::vector<ColourMatrix> Mvec(nsimd);

  vReal deltaS;
	vColourMatrix staple_sum1;  //sum of common terms
  // vColourMatrix staple_sum2;  //sum of common terms
  staple_sum1=zero;
  // staple_sum2=zero;

	for(int mu=0; mu<Nd; ++mu)
  {
    staple_sum1 += peekvCMatrix(U, x, mu) * adj(peekvg(g, shiftUp(g, x, mu)))
										+ adj(peekvg(g, shiftDown(g, x, mu)) *  peekvCMatrix(U, shiftDown(g, x, mu), mu));
    // staple_sum2 += peekvg(g, shiftDown(g, x, mu)) * peekvCMatrix(U, shiftDown(g, x, mu), mu);
  }

  for(int i=0; i<multi_hit; ++i)
  {
		vColourMatrix gx = g[x];
    for(int idx=0; idx<nsimd; ++idx)
      LieRandomize_serial(sRNG, Mvec[idx], scale);
    merge(vM, Mvec);

    gx_new = vM * gx;
		deltaS = - betaMM * toReal(TensorRemove( trace( (gx_new - gx) * staple_sum1 ) ));
                        // + toReal(TensorRemove( trace( adj(gx_new - gx) * staple_sum2 ) )) );
		// std::cout << "deltaS" << deltaS << std::endl;

    for(int idx=0; idx<nsimd; ++idx)
    {
      Real *pReal = (Real *) &deltaS;
      if(!metropolis_test(*(pReal+idx*2), sRNG)) Mvec[idx]=1.0;
    }

    merge(vM, Mvec);
    gx_new = vM * gx;
    g[x] = gx_new;
    // pokeIndex<LorentzIndex>(U[x], Uxmu_new, mu);
  }

}

void update_g_all(LatticeColourMatrix &g, const LatticeGaugeField &U, GridSerialRNG &sRNG, RealD betaMM, RealD scale, int multi_hit)
{
  for(int x=0; x<g.Grid()->oSites(); ++x)
  {
      update_gx(x, g, U, sRNG, betaMM, scale, multi_hit);
  }
}


Real GF_S(const LatticeColourMatrix &g, const LatticeGaugeField &Umu)
{
  LatticeColourMatrix s(Umu.Grid());
  s=zero;

  std::vector<LatticeColourMatrix> U(4, Umu.Grid());
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  for(int mu=0; mu<Nd; mu++)
    s += g * U[mu] * adj(Cshift(g, mu, 1));

  Complex p = TensorRemove( sum( trace(s) ));
  return - p.real();
}

void GF_metro(const LatticeGaugeField &U, LatticeColourMatrix &g, int nsweeps, Real betaMM, int multi_hit, double scale)
{
	int seed = std::chrono::system_clock::now().time_since_epoch().count();
	GridSerialRNG    sRNG;       sRNG.SeedFixedIntegers(std::vector<int>{seed});
	for(int i=0; i<nsweeps; ++i) update_g_all(g, U, sRNG, betaMM, scale, multi_hit);
}





}
