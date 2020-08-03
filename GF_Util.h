#include <cmath>



namespace Grid{

template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
  using namespace Grid;
  std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
  typename T::vector_object::scalar_object site;
  peekSite(site, field, coor);
  std::cout << site << std::endl;
}



template<class rtype,class ltype, int N>
strong_inline auto operator*(const iVector<rtype, N> &rhs, const iVector<ltype, N> &lhs) -> iVector<decltype(rhs(0)*lhs(0)), N>
{
	iVector<decltype(rhs(0)*lhs(0)), 4> ret;
	for(int c1=0;c1<N;c1++){
			mult(&ret._internal[c1],&rhs._internal[c1],&lhs._internal[c1]);
	}
	return ret;
}

namespace QCD{

using LatticeGaugeFieldSite = typename LatticeGaugeField::vector_object::scalar_object;

template <class T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &v) {
  if(v.empty()) {
    out << "[]";
    return out;
  }

  out << "[";
  for(size_t i=0; i<v.size()-1; ++i) out << v[i] << " "; // if v.size()=0, then v.size()-1 is a very large number (v.size() is unsigned!)
  out << v.back() << "]";
  return out;
}



void readField(LatticeGaugeField &U, const std::string &filename)
{
		FieldMetaData header;
    NerscIO::readConfiguration(U,header,filename);
}

void writeField(LatticeGaugeField &U, const std::string &filename)
{
		NerscIO::writeConfiguration(U,filename,0,0);
}


inline Coordinate operator-(const Coordinate &v1, const Coordinate &v2)
{
  int vsize = v1.size();
  Coordinate c(vsize);
  for(int i=0; i<vsize; ++i)
  {
    c[i] = v1[i] - v2[i];
  }
  return c;
}



inline std::vector<int> operator-(const std::vector<int> &v1, const std::vector<int> &v2)
{
  int vsize = v1.size();
  std::vector<int> c(vsize);
  for(int i=0; i<vsize; ++i)
  {
    c[i] = v1[i] - v2[i];
  }
  return c;
}

inline std::vector<int> operator+(const std::vector<int> &v1, const std::vector<int> &v2)
{
  int vsize = v1.size();
  std::vector<int> c(vsize);
  for(int i=0; i<vsize; ++i)
  {
    c[i] = v1[i] + v2[i];
  }
  return c;
}

inline std::vector<int> operator%(const std::vector<int> &v1, const std::vector<int> &v2)
{
  int vsize = v1.size();
  std::vector<int> c(vsize);
  for(int i=0; i<vsize; ++i)
  {
    c[i] = v1[i]%v2[i];
  }
  return c;
}



//return \sum_{x,\mu} Re\ tr[U_\mu(x)]
Real Omega_no_g(const LatticeGaugeField &Umu)
{
  Real ret;
  ret = WilsonLoops<PeriodicGimplR>::linkTrace(Umu) * Umu.Grid()->gSites() * 4.0 * 3.0;
  return ret;
}

//return \sum_{x,\mu} Re\ tr[g(x)U_\mu(x)g(x+\hat{\mu})]
Real Omega_g(const LatticeColourMatrix &g, const LatticeGaugeField &Umu)
{
  Real ret;
  LatticeGaugeField Utrans(Umu.Grid());
  Utrans = Umu;
  LatticeColourMatrix gg(g.Grid()); //in SU3::GaugeTransform, g is passed by non-const reference; FIXME
  gg = g;
  SU3::GaugeTransform(Utrans, gg);
  ret = WilsonLoops<PeriodicGimplR>::linkTrace(Utrans) * Utrans.Grid()->gSites() * 4.0 * 3.0;
  return ret;
}

//return Ta(U_\mu(x) g(x+\mu)^\dagger g(x))
LatticeGaugeField dOmegadU_g(const LatticeColourMatrix &g, const LatticeGaugeField &Umu)
{
  LatticeGaugeField ret(Umu.Grid());
  ret = Zero();

  LatticeColourMatrix s(Umu.Grid());

  std::vector<LatticeColourMatrix> U(4, Umu.Grid());
  for(int mu=0; mu<Nd; mu++) {
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);
		s = U[mu] * adj(Cshift(g, mu, 1)) * g;
    PokeIndex<LorentzIndex>(ret, s, mu);
  }

  ret = Ta(ret);
  return ret;
}

//return Ta(U_\mu(x) g(x+\mu)^\dagger g(x))
LatticeGaugeField dOmegadU_g(const LatticeColourMatrix &g, const std::vector<LatticeColourMatrix> &Umu)
{
  LatticeGaugeField ret(g.Grid());
  ret = Zero();

  LatticeColourMatrix s(g.Grid());

  for(int mu=0; mu<Nd; mu++) {
		s = Umu[mu] * adj(Cshift(g, mu, 1)) * g;
    PokeIndex<LorentzIndex>(ret, s, mu);
  }

  ret = Ta(ret);
  return ret;
}



Real dOmegaSquare2_no_g(const LatticeGaugeField &Umu)
{
  std::vector<LatticeColourMatrix> U(4, Umu.Grid());
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix s(Umu.Grid());
  s = Zero();
  for(int mu=0; mu<Nd; mu++)
    s += U[mu] - Cshift(U[mu], mu, -1);

  LatticeReal s1(Umu.Grid());
  LatticeReal s2(Umu.Grid());
  s2 = Zero();
  ColourMatrix ti;
  for(int i=0; i<8;++i)
  {
    SU<3>::generator(i, ti);
    s1 = toReal( timesI( trace(ti * s) ) );
    s2 += s1 * s1;
  }
  return TensorRemove(sum(s2));
}



Real dOmegaSquare2(const LatticeColourMatrix &g,const LatticeGaugeField &Umu)
{
  std::vector<LatticeColourMatrix> U(4, Umu.Grid());
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix s(Umu.Grid());
  s = Zero();
  for(int mu=0; mu<Nd; mu++)
    s += g * U[mu] * adj(Cshift(g, mu, 1)) - Cshift(g, mu, -1) * Cshift(U[mu], mu, -1) * adj(g);

  LatticeReal s1(Umu.Grid());
  LatticeReal s2(Umu.Grid());
  s2 = Zero();
  ColourMatrix ti;
  for(int i=0; i<8;++i)
  {
    SU<3>::generator(i, ti);
    s1 = toReal( timesI( trace(ti * s) ) );
    s2 += s1 * s1;
  }
  return TensorRemove(sum(s2));
}



// //assign a complex number to LatticeComplex
// void assign_lc(LatticeComplex &l, const std::complex<double> &s)
// {
// 	vComplex v;
// 	vComplex::conv_t con_v;
// 	for(int i=0; i<vComplex::Nsimd(); ++i) con_v.s[i] = s;
// 	v.v = con_v.v;
//
// 	parallel_for(int ss=0; ss<l.Grid()->oSites(); ++ss)
// 	{
// 		l._odata[ss] = v;
// 	}
// }
//
// //take eight component with respect to Gell-Mann matrices
// std::vector<Complex> peek_a(const ColourMatrix &m, const std::vector<ColourMatrix> &ta)
// {
//   std::vector<Complex> vm(8);
//   for(int a=0; a<8; ++a)  vm[a] = 2.0 * TensorRemove(trace(m * ta[a]));
//   return vm;
// }
//
// //For matrix at site x, LorentzIndex \mu, take eight component with respect to Gell-Mann matrices
// std::vector<Complex> peek_xa(const LatticeGaugeField &M, const std::vector<int> &coor, int mu,  const std::vector<ColourMatrix> &ta)
// {
//   iVector<iScalar<iMatrix<Complex, 3>>, 4> m;
//   peekSite(m, M, coor);
//   ColourMatrix m2;
//   m2() = m(mu); //cannot use peekIndex for non-vector datatype
//   return peek_a(m2, ta);
// }
//
//
// // mult iVector<xxx>; FIXME
// template<class rtype,class vtype,class mtype,int N>
// strong_inline void multVV(iVector<rtype,N> * __restrict__ ret,
//                  const iVector<vtype,N> * __restrict__ rhs,
//                  const iVector<mtype,N> * __restrict__ lhs){
//   for(int c1=0;c1<N;c1++){
//       mult(&ret->_internal[c1],&rhs->_internal[c1],&lhs->_internal[c1]);
//   }
// }


// //real(LatticeComplex) does not work; FIXME
// //LatticeComplex -> LatticeComplex
// template<class vobj> inline Lattice<vobj> zyd_real(const Lattice<vobj> &lhs){
//     Lattice<vobj> ret(lhs.Grid());
//     parallel_for(int ss=0;ss<lhs.Grid()->oSites();ss++){
//       ret._odata[ss] = real(lhs._odata[ss]);
//     }
//     return ret;
// };



// template <class T>
// void printRB(const Lattice<T> &lat) {
//   for(int i=0; i<lat._odata.size(); ++i) {
//     std::cout << lat[i] << std::endl;
//   }
// }

// double maxLattice(const LatticeReal& lat) {
//
//   typedef typename std::remove_reference<decltype(lat)>::type::scalar_type scalar_type;
//
//   double max_val = *((scalar_type *)& lat[0]()()());
//   #pragma omp parallel for reduction(max : max_val)
//   for(int ss=0;ss<lat.Grid()->oSites();ss++)
//   {
//     scalar_type *sobj = (scalar_type *)& lat[ss]()()();
//     for(int idx=0; idx<lat.Grid()->Nsimd(); ++idx)
//       if(*(sobj + idx) > max_val)
//           max_val = *(sobj + idx);
//   }
//
//   #ifndef GRID_COMMS_NONE
//   MPI_Allreduce(MPI_IN_PLACE, &max_val, 1, MPI_DOUBLE, MPI_MAX, lat.Grid()->communicator);
//   #endif
//
//   return max_val;
// }
//
// double minLattice(const LatticeReal& lat) {
//
//   typedef typename std::remove_reference<decltype(lat)>::type::scalar_type scalar_type;
//
//   double min_val = *((scalar_type *)& lat[0]()()());
//   #pragma omp parallel for reduction(min : min_val)
//   for(int ss=0;ss<lat.Grid()->oSites();ss++)
//   {
//     scalar_type *sobj = (scalar_type *)& lat[ss]()()();
//     for(int idx=0; idx<lat.Grid()->Nsimd(); ++idx)
//       if(*(sobj + idx) < min_val)
//           min_val = *(sobj + idx);
//   }
//
//   #ifndef GRID_COMMS_NONE
//   MPI_Allreduce(MPI_IN_PLACE, &min_val, 1, MPI_DOUBLE, MPI_MIN, lat.Grid()->communicator);
//   #endif
//
//   return min_val;
// }
//


}}
