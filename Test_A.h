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



void measure_A(const LatticeGaugeField &U, const std::vector<std::vector<int>> &coors) {

    // print_grid_field_site(U, {0,0,0,0});
    // print_grid_field_site(U, {1,0,0,0});
    LatticeGaugeField An(U._grid);
    An = timesMinusI(Log(U));
    // print_grid_field_site(An, {1,0,0,0});

    FFT theFFT((Grid::GridCartesian *)An._grid);
    theFFT.FFT_all_dim(An, An, FFT::forward);

    Lattice<iVector<iScalar<iVector<vComplex, 8 > >, 4> > alg(U._grid);

    for (int mu = 0; mu < Nd; mu++)
    {
      LatticeColourMatrix An_mu(An._grid);
      SU3::LatticeAlgebraVector alg_mu(U._grid);

      An_mu = peekLorentz(An, mu);
      SU3::projectOnAlgebra(alg_mu, An_mu);
      pokeLorentz(alg, alg_mu, mu);
    }
    // std::cout << alg << std::endl;

    for(const std::vector<int> &coor : coors) {
      // std::cout << "coor=" << coor << ": ";
      std::cout << "coor=";
      print_grid_field_site(alg, coor);
    }

}

}}
