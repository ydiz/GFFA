#pragma once

namespace Grid {

struct GFFAParams : Serializable {
public:
  bool isGFFA;
  bool isCell;
  double betaMM;

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(GFFAParams,
    bool, newHp,
    std::string, action,
    bool, add_fermion,
    double, beta,
    double, M,
    double, epsilon,   // When epsilon is set to 0., the force of zero mode is automatically set to 0; see GF_deltaU.h
    std::string, table_path,
    int, saveInterval,
    std::string, UFile,
    int, innerMC_N,
    int, hb_offset,
    std::vector<int>, cell_size,

    bool, measure_A,
    // double, fixed_P_k,
    std::vector<std::vector<int>>, measure_A_coors,

    double, maxHours // maximum allowed time in hours
  );

  template <class ReaderClass>
  GFFAParams(Reader<ReaderClass>& reader) {
    read(reader, "GFFA", *this);

    isGFFA = action.substr(0,2) == "GF";  // First two characters are "GF"
    isCell = action.substr(action.size() - 4) == "cell"; // Last 4 characters are "cell"
    betaMM = beta * M * M;

    std::cout << "isGFFA: " << isGFFA << std::endl;
    std::cout << "isCell: " << isCell << std::endl;
    // std::cout << *this << std::endl;
  }

};


struct My_DomainWallFermionParameters : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(My_DomainWallFermionParameters,
          RealD, mass,
          RealD, M5,
          unsigned int, Ls);

  template <class ReaderClass>
  My_DomainWallFermionParameters(Reader<ReaderClass>& reader) {
    read(reader, "DomainWallFermionParameters", *this);
    // std::cout << *this << std::endl;
  }
};

}

