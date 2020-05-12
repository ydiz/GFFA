#include <cmath>



namespace Grid{

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


// void printMem()
// {
// 	struct sysinfo myinfo;
// 	sysinfo(&myinfo); // there is a struct called sysinfo and a function called sysinfo as well
// 	double total_mem = myinfo.mem_unit * myinfo.totalram;
// 	total_mem /= (1024.*1024.);
// 	double free_mem = myinfo.mem_unit * myinfo.freeram;
// 	free_mem /= (1024.*1024.);
//
// 	printf("printMem node d: Memory: total: %.2f MB, avail: %.2f MB, used %.2f MB\n", total_mem, free_mem, total_mem-free_mem);
// }

void readField(LatticeGaugeField &U, const std::string &filename)
{
		FieldMetaData header;
			NerscIO::readConfiguration(U,header,filename);
}

void writeField(LatticeGaugeField &U, const std::string &filename)
{
		NerscIO::writeConfiguration(U,filename,0,0);
}


//assign a complex number to LatticeComplex
void assign_lc(LatticeComplex &l, const std::complex<double> &s)
{
	vComplex v;
	vComplex::conv_t con_v;
	for(int i=0; i<vComplex::Nsimd(); ++i) con_v.s[i] = s;
	v.v = con_v.v;

	parallel_for(int ss=0; ss<l._grid->oSites(); ++ss)
	{
		l._odata[ss] = v;
	}
}

//take eight component with respect to Gell-Mann matrices
std::vector<Complex> peek_a(const ColourMatrix &m, const std::vector<ColourMatrix> &ta)
{
  std::vector<Complex> vm(8);
  for(int a=0; a<8; ++a)  vm[a] = 2.0 * TensorRemove(trace(m * ta[a]));
  return vm;
}

//For matrix at site x, LorentzIndex \mu, take eight component with respect to Gell-Mann matrices
std::vector<Complex> peek_xa(const LatticeGaugeField &M, const std::vector<int> &coor, int mu,  const std::vector<ColourMatrix> &ta)
{
  iVector<iScalar<iMatrix<Complex, 3>>, 4> m;
  peekSite(m, M, coor);
  ColourMatrix m2;
  m2() = m(mu); //cannot use peekIndex for non-vector datatype
  return peek_a(m2, ta);
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


// mult iVector<xxx>; FIXME
template<class rtype,class vtype,class mtype,int N>
strong_inline void multVV(iVector<rtype,N> * __restrict__ ret,
                 const iVector<vtype,N> * __restrict__ rhs,
                 const iVector<mtype,N> * __restrict__ lhs){
  for(int c1=0;c1<N;c1++){
      mult(&ret->_internal[c1],&rhs->_internal[c1],&lhs->_internal[c1]);
  }
}


//real(LatticeComplex) does not work; FIXME
//LatticeComplex -> LatticeComplex
template<class vobj> inline Lattice<vobj> zyd_real(const Lattice<vobj> &lhs){
    Lattice<vobj> ret(lhs._grid);
    parallel_for(int ss=0;ss<lhs._grid->oSites();ss++){
      ret._odata[ss] = real(lhs._odata[ss]);
    }
    return ret;
};

//return \sum_{x,\mu} Re\ tr[U_\mu(x)]
Real Omega_no_g(const LatticeGaugeField &Umu)
{
  Real ret;
  ret = WilsonLoops<PeriodicGimplR>::linkTrace(Umu) * Umu._grid->gSites() * 4.0 * 3.0;
  return ret;
}

//return \sum_{x,\mu} Re\ tr[g(x)U_\mu(x)g(x+\hat{\mu})]
Real Omega_g(const LatticeColourMatrix &g, const LatticeGaugeField &Umu)
{
  Real ret;
  LatticeGaugeField Utrans(Umu._grid);
  Utrans = Umu;
  LatticeColourMatrix gg(g._grid); //in SU3::GaugeTransform, g is passed by non-const reference; FIXME
  gg = g;
  SU3::GaugeTransform(Utrans, gg);
  ret = WilsonLoops<PeriodicGimplR>::linkTrace(Utrans) * Utrans._grid->gSites() * 4.0 * 3.0;
  return ret;
}

//return Ta(U_\mu(x) g(x+\mu)^\dagger g(x))
LatticeGaugeField dOmegadU_g(const LatticeColourMatrix &g, const LatticeGaugeField &Umu)
{
  LatticeGaugeField ret(Umu._grid);
  ret = zero;

  LatticeColourMatrix s(Umu._grid);

  std::vector<LatticeColourMatrix> U(4, Umu._grid);
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
  LatticeGaugeField ret(g._grid);
  ret = zero;

  LatticeColourMatrix s(g._grid);

  for(int mu=0; mu<Nd; mu++) {
		s = Umu[mu] * adj(Cshift(g, mu, 1)) * g;
    PokeIndex<LorentzIndex>(ret, s, mu);
  }

  ret = Ta(ret);
  return ret;
}



Real dOmegaSquare2_no_g(const LatticeGaugeField &Umu)
{
  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix s(Umu._grid);
  s = zero;
  for(int mu=0; mu<Nd; mu++)
    s += U[mu] - Cshift(U[mu], mu, -1);

  LatticeReal s1(Umu._grid);
  LatticeReal s2(Umu._grid);
  s2 = zero;
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
  std::vector<LatticeColourMatrix> U(4, Umu._grid);
  for(int mu=0; mu<Nd; mu++)
    U[mu] = PeekIndex<LorentzIndex>(Umu, mu);

  LatticeColourMatrix s(Umu._grid);
  s = zero;
  for(int mu=0; mu<Nd; mu++)
    s += g * U[mu] * adj(Cshift(g, mu, 1)) - Cshift(g, mu, -1) * Cshift(U[mu], mu, -1) * adj(g);

  LatticeReal s1(Umu._grid);
  LatticeReal s2(Umu._grid);
  s2 = zero;
  ColourMatrix ti;
  for(int i=0; i<8;++i)
  {
    SU<3>::generator(i, ti);
    s1 = toReal( timesI( trace(ti * s) ) );
    s2 += s1 * s1;
  }
  return TensorRemove(sum(s2));
}

template <class T>
void printRB(const Lattice<T> &lat) {
  for(int i=0; i<lat._odata.size(); ++i) {
    std::cout << lat[i] << std::endl;
  }
}

double maxLattice(const LatticeReal& lat) {

  typedef typename std::remove_reference<decltype(lat)>::type::scalar_type scalar_type;

  double max_val = *((scalar_type *)& lat[0]()()());
  #pragma omp parallel for reduction(max : max_val)
  for(int ss=0;ss<lat._grid->oSites();ss++)
  {
    scalar_type *sobj = (scalar_type *)& lat[ss]()()();
    for(int idx=0; idx<lat._grid->Nsimd(); ++idx)
      if(*(sobj + idx) > max_val)
          max_val = *(sobj + idx);
  }

  #ifndef GRID_COMMS_NONE
  MPI_Allreduce(MPI_IN_PLACE, &max_val, 1, MPI_DOUBLE, MPI_MAX, lat._grid->communicator);
  #endif

  return max_val;
}

double minLattice(const LatticeReal& lat) {

  typedef typename std::remove_reference<decltype(lat)>::type::scalar_type scalar_type;

  double min_val = *((scalar_type *)& lat[0]()()());
  #pragma omp parallel for reduction(min : min_val)
  for(int ss=0;ss<lat._grid->oSites();ss++)
  {
    scalar_type *sobj = (scalar_type *)& lat[ss]()()();
    for(int idx=0; idx<lat._grid->Nsimd(); ++idx)
      if(*(sobj + idx) < min_val)
          min_val = *(sobj + idx);
  }

  #ifndef GRID_COMMS_NONE
  MPI_Allreduce(MPI_IN_PLACE, &min_val, 1, MPI_DOUBLE, MPI_MIN, lat._grid->communicator);
  #endif

  return min_val;
}



}}
