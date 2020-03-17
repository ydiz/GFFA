#include <Grid/Grid.h>

#include "Test_gauge_force.h"

namespace Grid {
namespace QCD {


template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
strong_inline iMatrix<vtype,N> Log(const iMatrix<vtype,N> &arg, Integer Nlog = 12 )
{
  typedef iMatrix<vtype,N> mat;
  mat unit(1.0);
  mat argMinusOne = arg - unit;

  mat argMinusOnePowerI = unit;
  mat temp = zero;
  for(int i=1; i<=Nlog; ++i){
    argMinusOnePowerI *= argMinusOne;
    if(i%2==1) temp += argMinusOnePowerI * (1. / double(i));
    else temp -= argMinusOnePowerI * (1. / double(i));
  }
  return temp;

}


LatticeGaugeField Log(const LatticeGaugeField &lat) {
  LatticeGaugeField rst(lat._grid);

  parallel_for(int ss=0; ss<lat._grid->oSites(); ss++) {
    for(int mu=0; mu<4; ++mu) {
      rst[ss](mu)() = Log(lat[ss](mu)());
    }
  }

  return rst;
}



void measure_A(const LatticeGaugeField &U, const std::vector<std::vector<int>> &coors, bool log=true) {

  LatticeGaugeField An(U._grid);
  if(log) An = timesMinusI(Log(U));
  else An = U;

  FFT theFFT((Grid::GridCartesian *)An._grid);
  // std::cout << "measure_A" << std::endl;
  // print_grid_field_site(An, {0,0,0,0});
  // print_grid_field_site(An, {1,1,1,1});
  theFFT.FFT_all_dim(An, An, FFT::forward);
  // print_grid_field_site(An, {0,0,0,0});
  // print_grid_field_site(An, {1,1,1,1});
  double vol = 1.0;
  for(int d=0; d<4; ++d)   vol = vol * U._grid->_fdimensions[d];
  An = An * (1. / std::sqrt(vol));

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


}}
