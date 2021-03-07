#pragma once

namespace Grid {

using vLorentzScalar = iVector<iScalar<iScalar<vComplex> >, Nd >;
using LatticeLorentzScalar = Lattice<vLorentzScalar>;

LatticeLorentzScalar get_mask(GridBase *grid, const Coordinate &cell_size) { // grid for entire lattice, not for inner cell

  LatticeLorentzScalar mask(grid); mask = 1.;
  LatticeLorentzScalar zero(grid); zero = 0.;
  for(int mu=0; mu<4; mu++){
    LatticeInteger coor_mu(grid);
    LatticeCoordinate(coor_mu, mu);

    mask = where(coor_mu >= (Integer)cell_size[mu], zero, mask);   // Must explicitly convert to Integer to compile.
  }
  return mask;
}

}


