#pragma once

namespace Grid {

using vLorentzScalar = iVector<iScalar<iScalar<vComplex> >, Nd >;
using LatticeLorentzScalar = Lattice<vLorentzScalar>;

// set a mask on every node
LatticeLorentzScalar get_mask(GridBase *grid, const Coordinate &cell_size) { // grid for entire lattice, not for inner cell

  LatticeLorentzScalar mask(grid);

  Coordinate ldf = grid->_ldimensions;

  autoView(mask_v, mask, CpuWrite);
  // thread_for(ss, grid->lSites(), {
  accelerator_for(ss, grid->lSites(), 1, {
    Coordinate lcoor;
    // Coordinate lcoor, gcoor;
    // localIndexToLocalGlobalCoor(grid, ss, lcoor, gcoor);
    Lexicographic::CoorFromIndex(lcoor, ss, ldf);

    typename LatticeLorentzScalar::vector_object::scalar_object m;
    m = 1.0;
    for(int mu=0; mu<4; ++mu) {
      if(lcoor[mu] >= cell_size[mu]) {  // If x is outside cell, m(mu) = 0, for any mu
        m = 0.;
        break;
      }
      if(lcoor[mu] == cell_size[mu] - 1) m(mu)()() = 0.; // [1,5,3,5] ->  mask = (1,0,1,0) // If x is outside cell, m(mu) = 0, for any mu
    }

    pokeLocalSite(m, mask_v, lcoor);
  });
  return mask;
}


// FIXME: set a mask on every node
LatticeLorentzScalar get_cell_mask(GridBase *cell_grid) { // grid for entire lattice, not for inner cell

  LatticeLorentzScalar mask(cell_grid);

  Coordinate ldf = cell_grid->_ldimensions;

  autoView(mask_v, mask, CpuWrite);
  thread_for(ss, cell_grid->lSites(), {
    // Coordinate lcoor, gcoor;
    // localIndexToLocalGlobalCoor(cell_grid, ss, lcoor, gcoor);
    Coordinate lcoor;
    // Coordinate lcoor, gcoor;
    // localIndexToLocalGlobalCoor(grid, ss, lcoor, gcoor);
    Lexicographic::CoorFromIndex(lcoor, ss, ldf);

    typename LatticeLorentzScalar::vector_object::scalar_object m;
    m = 1.0;
    for(int mu=0; mu<4; ++mu) {
      if(lcoor[mu] == ldf[mu] - 1) m(mu)()() = 0.; // [1,5,3,5] ->  mask = (1,0,1,0) 
    }
    pokeLocalSite(m, mask_v, lcoor);
  });
  return mask;
}

// LatticeLorentzScalar get_mask(GridBase *grid, const Coordinate &cell_size) { // grid for entire lattice, not for inner cell
//
//   LatticeLorentzScalar mask(grid);
//
//   autoView(mask_v, mask, CpuWrite);
//   thread_for(ss, grid->lSites(), {
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(grid, ss, lcoor, gcoor);
//
//     typename LatticeLorentzScalar::vector_object::scalar_object m;
//     m = 1.0;
//     for(int mu=0; mu<4; ++mu) {
//       if(gcoor[mu] >= cell_size[mu]) {  // If x is outside cell, m(mu) = 0, for any mu
//         m = 0.;
//         break;
//       }
//       if(gcoor[mu] == cell_size[mu] - 1) m(mu)()() = 0.; // [1,5,3,5] ->  mask = (1,0,1,0) // If x is outside cell, m(mu) = 0, for any mu
//     }
//
//     pokeLocalSite(m, mask_v, lcoor);
//   });
//   return mask;
// }
//
//
// LatticeLorentzScalar get_cell_mask(GridBase *cell_grid) { // grid for entire lattice, not for inner cell
//
//   LatticeLorentzScalar mask(cell_grid);
//
//   Coordinate cell_size = cell_grid->_fdimensions;
//
//   autoView(mask_v, mask, CpuWrite);
//   thread_for(ss, cell_grid->lSites(), {
//     Coordinate lcoor, gcoor;
//     localIndexToLocalGlobalCoor(cell_grid, ss, lcoor, gcoor);
//
//     typename LatticeLorentzScalar::vector_object::scalar_object m;
//     m = 1.0;
//     for(int mu=0; mu<4; ++mu) {
//       if(gcoor[mu] == cell_size[mu] - 1) m(mu)()() = 0.; // [1,5,3,5] ->  mask = (1,0,1,0) 
//     }
//     pokeLocalSite(m, mask_v, lcoor);
//   });
//   return mask;
// }



}


