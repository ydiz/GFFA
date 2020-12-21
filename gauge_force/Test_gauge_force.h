#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

int smod(const int x, const int len)
{
  if (x * 2 <= len) return x;
  else return x - len;
}

// std::vector<int> smod(const std::vector<int> &v, const std::vector<int> &lens)
Coordinate smod(const Coordinate &v, const Coordinate &lens)
{
  Coordinate rst(v.size());
  for(size_t i=0; i<v.size(); ++i) rst[i] = smod(v[i], lens[i]);
  return rst;  
}



// get (|k|, average F(|k|))
void aggregate(const LatticeGaugeField &force, std::vector<double> &sums, std::vector<double> &counts, double interval = 1.) {
  Coordinate fdims = force.Grid()->_fdimensions;
  double max_r = std::sqrt(0.25 * (fdims[0]*fdims[0] + fdims[1]*fdims[1] + fdims[2]*fdims[2] + fdims[3]*fdims[3]));

  sums.resize(int(max_r / interval)+1);
  counts.resize(int(max_r / interval)+1);
  // std::cout << "rst size: " << rst.size() << std::endl;
  

  autoView(force_v, force, CpuRead);
  
  // FIXME: this is wrong!!!! Will be data race in parallel for; need to use "reduction" in "parallel for"
  // https://stackoverflow.com/questions/43168661/openmp-and-reduction-on-stdvector
  // parallel_for(int ss=0; ss<force._grid->lSites(); ss++) {
  for(int ss=0; ss<force.Grid()->lSites(); ss++) {  // FIXME: I disabled parallelism
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(force.Grid(), ss, lcoor, gcoor);

    gcoor = smod(gcoor, fdims);
    double r = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2] + gcoor[3]*gcoor[3]);

    typename LatticeGaugeField::vector_object::scalar_object m;
    peekLocalSite(m, force_v, lcoor);

    // double val = norm2(real(m));  // FIXME: Dec 05, 2020: why take real part?
    double val = norm2(m); // Dec 10.
    // assert(int(r)<=16); 
    sums[int(r / interval)] += val;
    counts[int(r / interval)] += 1.0;


    if( abs(r - 1.) < 0.0001 ) {
      std::cout << "gcoor: " << gcoor << std::endl;
      std::cout << "sum: " << sums[10] << std::endl;
      std::cout << "count: " << counts[10] << std::endl;
    }

  }

  force.Grid()->GlobalSumVector(sums.data(), sums.size()); // sum over nodes
  force.Grid()->GlobalSumVector(counts.data(), counts.size());
}

std::vector<double> get_average_force(const LatticeGaugeField &force, double interval) {
  std::vector<double> sums;
  std::vector<double> counts;
  aggregate(force, sums, counts, interval);
  std::cout << "counts: " << counts << std::endl;

  // std::cout << sums << std::endl;
  // std::cout << counts << std::endl;

  std::vector<double> avg(sums.size());
  for(int i=0; i<sums.size(); ++i) avg[i] = (counts[i]!=0 ? sums[i] / counts[i] : 0.);

  return avg;
}


LatticeGaugeField PL_projection(const LatticeGaugeField &P, double epsilon) {
  LatticeGaugeField PL(P.Grid());

  static Momenta_k KK(P.Grid(), 0., epsilon, true);
  
  LatticeColourMatrix sinKExpDotPk(P.Grid());
  // sinKExpDotPk = KK.sinKPsExpDotP_func(P);
  sinKExpDotPk = KK.sinKNgExpDotP_func(P);

  LatticeColourMatrix PLmu(P.Grid());
  for(int mu=0; mu<Nd; ++mu) {
    // PLmu = KK.sinKNgExp[mu] * sinKExpDotPk;
    PLmu = KK.sinKPsExp[mu] * sinKExpDotPk;
    pokeLorentz(PL, PLmu, mu);
  }

  // PL = KK.one / KK.sinKEpsilonSquare * PL;
  PL = KK.one / KK.FourSinKSquareEpsilon * PL;

  return PL;
}

std::ostream& operator<<(std::ostream &out, const std::vector<double> &vec) {
  out << "[";
  for(auto x: vec) out << x << " ";
  out << "]";
  return out;
}



void get_force_stats(LatticeGaugeField &force, double interval, double epsilon) {
  std::cout << "before Fourier transformation" << std::endl;
  print_grid_field_site(force, {1,0,0,0});
  print_grid_field_site(force, {0,1,0,0});


  Coordinate fdims = force.Grid()->_fdimensions;
  double V = fdims[0] * fdims[1] * fdims[2] * fdims[3];

  FFT theFFT((Grid::GridCartesian *)force.Grid());

  // LatticeGaugeField force_k(force.Grid());   // does not make a difference
  // theFFT.FFT_all_dim(force_k, force, FFT::forward);
  // force_k *= 1. / sqrt(V);
  // std::cout << "after Fourier transformation" << std::endl;
  // print_grid_field_site(force_k, {1,0,0,0});
  // print_grid_field_site(force_k, {0,1,0,0});

  theFFT.FFT_all_dim(force, force, FFT::forward);
  force *= 1. / sqrt(V);

  std::cout << "after Fourier transformation" << std::endl;
  print_grid_field_site(force, {1,0,0,0});
  print_grid_field_site(force, {0,1,0,0});

  LatticeGaugeField force_L(force.Grid());
  LatticeGaugeField force_T(force.Grid());
  force_L = PL_projection(force, epsilon);
  force_T = force - force_L;

  std::cout << "force avg: " << get_average_force(force, interval) << std::endl;
  // std::cout << "force_L avg: " << get_average_force(force_L, interval) << std::endl;  // FIXME: uncomment this
  // std::cout << "force_T avg: " << get_average_force(force_T, interval) << std::endl;
}

}}
