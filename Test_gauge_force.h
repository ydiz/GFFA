#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

void localIndexToLocalGlobalCoor(GridBase *grid, int ss, std::vector<int> &lcoor, std::vector<int> &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  std::vector<int> processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}

int smod(const int x, const int len)
{
  if (x * 2 <= len) return x;
  else return x - len;
}

std::vector<int> smod(const std::vector<int> &v, const std::vector<int> &lens)
{
  std::vector<int> rst(v.size());
  for(size_t i=0; i<v.size(); ++i) rst[i] = smod(v[i], lens[i]);
  return rst;  
}



// get (|k|, average F(|k|))
void aggregate(const LatticeGaugeField &force, std::vector<double> &sums, std::vector<double> &counts, double interval = 1.) {
  std::vector<int> fdims = force._grid->_fdimensions;
  double max_r = std::sqrt(0.25 * (fdims[0]*fdims[0] + fdims[1]*fdims[1] + fdims[2]*fdims[2] + fdims[3]*fdims[3]));

  sums.resize(int(max_r / interval)+1);
  counts.resize(int(max_r / interval)+1);
  // std::cout << "rst size: " << rst.size() << std::endl;
  parallel_for(int ss=0; ss<force._grid->lSites(); ss++) {
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(force._grid, ss, lcoor, gcoor);

    gcoor = smod(gcoor, fdims);
    double r = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2] + gcoor[3]*gcoor[3]);

    typename LatticeGaugeField::vector_object::scalar_object m;
    peekLocalSite(m, force, lcoor);

    double val = norm2(real(m));
    // double val = norm2(m);
    // assert(int(r)<=16);
    sums[int(r / interval)] += val;
    counts[int(r / interval)] += 1.0;
  }
}

std::vector<double> get_average_force(const LatticeGaugeField &force, double interval) {
  std::vector<double> sums;
  std::vector<double> counts;
  aggregate(force, sums, counts, interval);

  force._grid->GlobalSumVector(sums.data(), sums.size());
  force._grid->GlobalSumVector(counts.data(), counts.size());
  // std::cout << sums << std::endl;
  // std::cout << counts << std::endl;

  std::vector<double> avg(sums.size());
  for(int i=0; i<sums.size(); ++i) avg[i] = (counts[i]!=0 ? sums[i] / counts[i] : 0.);

  return avg;
}


LatticeGaugeField PL_projection(const LatticeGaugeField &P, double epsilon) {
  LatticeGaugeField PL(P._grid);

  static Momenta_k KK(P._grid, 0., epsilon, true);
  
  LatticeColourMatrix sinKExpDotPk(P._grid);
  sinKExpDotPk = KK.sinKPsExpDotP_func(P);

  LatticeColourMatrix PLmu(P._grid);
  for(int mu=0; mu<Nd; ++mu) {
    PLmu = KK.sinKNgExp[mu] * sinKExpDotPk;
    pokeLorentz(PL, PLmu, mu);
  }

  PL = KK.one / KK.sinKEpsilonSquare * PL;

  return PL;
}

std::ostream& operator<<(std::ostream &out, const std::vector<double> &vec) {
  out << "[";
  for(auto x: vec) out << x << " ";
  out << "]";
  return out;
}

template<class T>
void print_grid_field_site(const T &field, const std::vector<int> &coor) {
  using namespace Grid;
  std::cout << coor << std::endl;
  // std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
  typename T::vector_object::scalar_object site;
  peekSite(site, field, coor);
  std::cout << site << std::endl;
}


void get_force_stats(LatticeGaugeField &force, double interval, double epsilon) {
  FFT theFFT((Grid::GridCartesian *)force._grid);
  theFFT.FFT_all_dim(force, force, FFT::forward);

  // print_grid_field_site(force, {0,0,0,0});
  
  LatticeGaugeField force_L(force._grid);
  LatticeGaugeField force_T(force._grid);
  force_L = PL_projection(force, epsilon);
  force_T = force - force_L;

  std::cout << "force avg: " << get_average_force(force, interval) << std::endl;
  std::cout << "force_L avg: " << get_average_force(force_L, interval) << std::endl;
  std::cout << "force_T avg: " << get_average_force(force_T, interval) << std::endl;
}

}}
