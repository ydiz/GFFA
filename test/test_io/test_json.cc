#include <Grid/Grid.h>


using namespace std;
using namespace Grid;

struct GFFAParams : Serializable {
public:
  std::vector<int> fdims; // lattice size

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GFFAParams,
    int, traj,
    std::vector<int>, coor,
    std::vector<std::vector<int>>, coor2
    // std::string, coor
  );

  template <class ReaderClass>
  GFFAParams(Reader<ReaderClass>& reader) {
    read(reader, "GFFA", *this);
    // GridCmdOptionIntVector(lat, fdims);
    coor2[0][0] = 10;
    std::cout << coor2 << std::endl;
    std::cout << *this << std::endl;
  }

};





int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  JSONReader reader("GFFA.json");
  GFFAParams GFFA_params(reader);


  // HMCparameters hmc_params;
  //
  // JSONWriter writer("hmc.json");
  // write(writer, "testObject", hmc_params);
  // std::cout << "after writing" << std::endl;

  Grid_finalize();

} // main
