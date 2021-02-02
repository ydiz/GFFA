#pragma once

namespace Grid {

// arg is between [-pi, pi]
void lat_cmpl_to_arg_norm(const LatticeComplex &lat, LatticeComplex &arg, LatticeComplex &norm) {

  // pokeLocalSite does not work for LatticeReal; Assertion `sizeof(sobj)*Nsimd == sizeof(vobj)' failed.

  autoView(lat_v, lat, CpuRead);
  autoView(arg_v, arg, CpuWrite);
  autoView(norm_v, norm, CpuWrite);

  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor;
    lat.Grid()->LocalIndexToLocalCoor(ss, lcoor);

    typename LatticeComplex::vector_object::scalar_object m, m_arg, m_norm;
    peekLocalSite(m, lat_v, lcoor);

    m_arg()()() = Complex(std::arg(m()()()), 0.);
    pokeLocalSite(m_arg, arg_v, lcoor);

    m_norm()()() = Complex(std::abs(m()()()), 0.);  // abs(complex) returns its magnitude; 
    pokeLocalSite(m_norm, norm_v, lcoor);
  });
}

void avg_std_phases_norms(const LatticeComplex &lat, double &avg, double &std) {
  Coordinate fdims = lat.Grid()->FullDimensions();
  int V = 1;
  for(int mu=0; mu<4; ++mu) V *= fdims[mu];

  LatticeReal lat_r = toReal(lat);

  avg = sum(lat_r)()()() / double(V);

  LatticeReal tmp = lat_r - avg;
  LatticeReal tmp2 = tmp * tmp;
  std = std::sqrt(sum(tmp2)()()() / double(V));
}

std::vector<uint64_t> binning_phases(const LatticeComplex &lat) {
  std::vector<uint64_t> counts(64); // [-3.2, 3.2); interval is 0.1 // GlobalSumVector only works for uint64_t, not int

  autoView(lat_v, lat, CpuRead);
  //  Cannot use thread parallelism; will have data race
  for(int ss=0; ss<lat.Grid()->lSites(); ss++) {  
    Coordinate lcoor;
    lat.Grid()->LocalIndexToLocalCoor(ss, lcoor);

    typename LatticeComplex::vector_object::scalar_object m;
    peekLocalSite(m, lat_v, lcoor);
    double phase = m()()().real();

    int idx = int(phase + 3.2 / 0.1);  // first bin is [-3.2  -3.1)
    counts.at(idx) += 1;
  }
  lat.Grid()->GlobalSumVector(counts.data(), counts.size());
  return counts;
}

std::vector<uint64_t> binning_norms(const LatticeComplex &lat) {
  std::vector<uint64_t> counts(31); // [0, 3]; interval is 0.1 // GlobalSumVector only works for uint64_t, not int

  autoView(lat_v, lat, CpuRead);
  //  Cannot use thread parallelism; will have data race
  for(int ss=0; ss<lat.Grid()->lSites(); ss++) {  
    Coordinate lcoor;
    lat.Grid()->LocalIndexToLocalCoor(ss, lcoor);

    typename LatticeComplex::vector_object::scalar_object m;
    peekLocalSite(m, lat_v, lcoor);
    double norm = m()()().real();

    int idx = int(norm / 0.1);  // first bin is [0 0.1)
    counts.at(idx) += 1;
  }
  lat.Grid()->GlobalSumVector(counts.data(), counts.size());
  return counts;
}




}
