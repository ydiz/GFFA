#include <Grid/Grid.h>


namespace Grid{


template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
  using namespace Grid;
  std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
  typename T::vector_object::scalar_object site;
  peekSite(site, field, coor);
  std::cout << site << std::endl;
}

using LatticeGaugeFieldSite = typename LatticeGaugeField::vector_object::scalar_object;

using Field = LatticeGaugeField;

void localIndexToLocalGlobalCoor(GridBase *grid, int ss, Coordinate &lcoor, Coordinate &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  Coordinate processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}


void TrajectoryComplete(int traj, Field &U) {

  Coordinate fdims = U.Grid()->FullDimensions();
  int V = fdims[0] * fdims[1] * fdims[2];
  int T = fdims[3];

  // for(int mu=0; mu<4; ++mu) {       // Direction of Polaykov loop
  for(int mu=3; mu<4; ++mu) {       // Direction of Polaykov loop
    LatticeColourMatrix Umu = peekLorentz(U, mu); // U_mu
    LatticeColourMatrix Umu_shift = Umu;
    LatticeColourMatrix tmp = Umu;
   
    // for(int i=0; i<fdims[mu]-1; ++i) tmp = tmp * Cshift(tmp, mu, 1);
    for(int i=0; i<fdims[mu]-1; ++i) {
      Umu_shift = Cshift(Umu_shift, mu, 1);
      tmp = tmp * Umu_shift;
    }

    // get three directions perpendicular to direction mu
    std::vector<int> perp_dirs {0,1,2,3};
    perp_dirs.erase(perp_dirs.begin() + mu);   
    int dir1 = perp_dirs[0], dir2 = perp_dirs[1], dir3 = perp_dirs[2];

    std::vector<Complex> polya(V);

    autoView(tmp_v , tmp, CpuRead);
    thread_for(ss, U.Grid()->lSites(), {
      Coordinate lcoor, gcoor;
      localIndexToLocalGlobalCoor(U.Grid(), ss, lcoor, gcoor);

      if(gcoor[mu] == 0) {
        typename LatticeColourMatrix::vector_object::scalar_object m;
        peekLocalSite(m, tmp_v, lcoor);

        typename LatticeComplex::vector_object::scalar_object m2;
        m2 = trace(m);

        int idx = gcoor[dir1] * fdims[dir1] * fdims[dir2] + gcoor[dir2] * fdims[dir2] + gcoor[dir3];
        polya[idx] = m2()()();
      }
    });
    int def_prec = std::cout.precision();
    std::cout << std::setprecision(3) << "Polyakov loop: [ " << traj << " ] mu = " << mu << ": " << polya << std::endl;

    // Complex avg_polya = 0.;
    // for(int i=0; i<polya.size(); ++i) avg_polya += polya[i];
    // avg_polya /= double(polya.size());
    //
    // std::cout << GridLogMessage
    //     << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
    //     << "AvgPolyakovLoop: [ " << traj << " ] mu = " << mu << ": " << avg_polya << std::endl;
    // std::cout.precision(def_prec);
  }

}




}

using namespace std;
using namespace Grid;

int main(int argc, char** argv) {
	Grid_init(&argc, &argv);

	Coordinate simd_layout = GridDefaultSimd(4,vComplex::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
	Coordinate latt_size  ({8,8,8,8});

	GridCartesian Grid(latt_size, simd_layout, mpi_layout);
	LatticeGaugeField U(&Grid);
	FieldMetaData header;

	NerscIO::readConfiguration(U, header, "./ckpoint_lat.1000");
  std::cout << "U site {0,0,0,0}" << std::endl;
  print_grid_field_site(U, {0,0,0,0});


  iMatrix<Complex, 3> rst = Zero();
  rst(0,0) = 1.; 
  rst(1,1) = 1.; 
  rst(2,2) = 1.; 
  for(int t=0; t<8; ++t) {
    LatticeGaugeFieldSite m;
    peekSite(m, U, Coordinate({0,0,0,t}));
    std::cout << m(3)() << std::endl;
    rst = rst * m(3)();
  }
  std::cout << "rst: " << rst << std::endl;
  std::cout << "trace: " << trace(rst) << std::endl;

  TrajectoryComplete(1000, U);




	// cout << U << endl;	

	return 0;
}
