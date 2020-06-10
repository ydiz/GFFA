#include <Grid/Grid.h>

namespace Grid {
namespace QCD {



std::vector<double> operator+(const std::vector<double> &x1, const std::vector<double> &x2) {
  std::vector<double> ret(x1.size());
  for(int i=0; i<ret.size(); ++i) ret[i] = x1[i] + x2[i];
  return ret;
}

std::vector<double> operator*(double x, const std::vector<double> &vec) {
  std::vector<double> ret(vec.size());
  for(int i=0; i<ret.size(); ++i) ret[i] = x * vec[i];
  return ret;
}


double maxNorm(const LatticeGaugeField& U) {
  Lattice<iVector<iScalar<iScalar<vRealD>>, Nd> > U_norm2(U._grid);

  parallel_for(int ss=0;ss<U_norm2._grid->oSites();ss++){
    for(int mu=0; mu<Nd; mu++) {
      U_norm2[ss](mu)()() = 0.;
      for(int c1=0;c1<Nc;c1++)
        for(int c2=0;c2<Nc;c2++)
            // U_norm2[ss](mu)()() += toReal(U[ss](mu)()(c1, c2) * U[ss](mu)()(c1, c2));
            U_norm2[ss](mu)()() += toReal(conjugate(U[ss](mu)()(c1, c2)) * U[ss](mu)()(c1, c2));
    }
  }

  typedef typename decltype(U_norm2)::scalar_type scalar_type;

  double max_val = 0.;
  #pragma omp parallel for reduction(max : max_val)
  for(int ss=0;ss<U_norm2._grid->oSites();ss++)
  {
    for(int mu=0; mu<Nd; mu++) {
      scalar_type *sobj = (scalar_type *)& U_norm2[ss](mu)()();
      for(int idx=0; idx<U_norm2._grid->Nsimd(); ++idx)
        if( *(sobj + idx) > max_val)
            max_val = *(sobj + idx);
    }
  }

  #ifndef GRID_COMMS_NONE
  MPI_Allreduce(MPI_IN_PLACE, &max_val, 1, MPI_DOUBLE, MPI_MAX, U_norm2._grid->communicator);
  #endif

  return max_val;
}


template<class T>
void writeVector(const std::vector<T> &vec, int traj, double tau, const std::string &filename, int rank) {
  if(rank == 0) {
    std::ofstream outFile(filename, std::ofstream::app);
    outFile << "traj " << traj << " tau " << tau << ": ";
    for(const T& x : vec) outFile << x << " ";
    outFile << "\n";
    outFile.close();
  }
#ifdef GRID_COMMS_MPI3
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


// void readField(LatticeGaugeField &U, const std::string &filename)
// {
// 		FieldMetaData header;
// 			NerscIO::readConfiguration(U,header,filename);
// }
//
// void writeField(LatticeGaugeField &U, const std::string &filename)
// {
// 		NerscIO::writeConfiguration(U,filename,0,0);
// }
//


}
}
