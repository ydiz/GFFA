namespace Grid{
namespace QCD{

template <class vcplx>
strong_inline vcplx matrix_mult_elem(const SU3::iSUnMatrix<vcplx> &m, const SU3::iSUnMatrix<vcplx> &n, int i0, int i1) {
  return m()()(i0, 0) * n()()(0, i1) + m()()(i0, 1) * n()()(1, i1) + m()()(i0, 2) * n()()(2, i1);
}

template <class vcplx>
void my_su2Extract(Lattice<iSinglet<vcplx> > &Determinant,
                       Lattice<SU3::iSU2Matrix<vcplx> > &subgroup,
                       const Lattice<SU3::iSUnMatrix<vcplx> > &link,
                       const Lattice<SU3::iSUnMatrix<vcplx> > &staple,
                       int su2_index) {
  GridBase *grid(link._grid);

  static int subgroup_index[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  int i0 = subgroup_index[su2_index][0];
  int i1 = subgroup_index[su2_index][1];

  parallel_for (int ss = 0; ss < grid->oSites(); ss++) {
    subgroup._odata[ss]()()(0, 0) = matrix_mult_elem(link[ss], staple[ss], i0, i0);
    subgroup._odata[ss]()()(0, 1) = matrix_mult_elem(link[ss], staple[ss], i0, i1);
    subgroup._odata[ss]()()(1, 0) = matrix_mult_elem(link[ss], staple[ss], i1, i0);
    subgroup._odata[ss]()()(1, 1) = matrix_mult_elem(link[ss], staple[ss], i1, i1);

    subgroup._odata[ss] = 0.5 * (subgroup._odata[ss] - adj(subgroup._odata[ss]) + trace(adj(subgroup._odata[ss])));

    Determinant._odata[ss] = subgroup._odata[ss]()()(0, 0) * subgroup._odata[ss]()()(1, 1)
                            - subgroup._odata[ss]()()(0, 1) * subgroup._odata[ss]()()(1, 0);
  }
}

template <class vcplx>
void my_su2Insert(const Lattice<SU3::iSU2Matrix<vcplx> > &subgroup,
                      Lattice<SU3::iSUnMatrix<vcplx> > &link, int su2_index) {
  GridBase *grid(link._grid);

  static int subgroup_index[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  int i0 = subgroup_index[su2_index][0];
  int i1 = subgroup_index[su2_index][1];

  parallel_for (int ss = 0; ss < grid->oSites(); ss++) {
    vcplx link_elem[2][3];
    for(int i=0; i<3; ++i) {link_elem[0][i] = link._odata[ss]()()(i0, i); link_elem[1][i] = link._odata[ss]()()(i1, i);}

    link._odata[ss]()()(i0, 0) = subgroup._odata[ss]()()(0, 0) * link_elem[0][0] + subgroup._odata[ss]()()(0, 1) * link_elem[1][0];
    link._odata[ss]()()(i0, 1) = subgroup._odata[ss]()()(0, 0) * link_elem[0][1] + subgroup._odata[ss]()()(0, 1) * link_elem[1][1];
    link._odata[ss]()()(i0, 2) = subgroup._odata[ss]()()(0, 0) * link_elem[0][2] + subgroup._odata[ss]()()(0, 1) * link_elem[1][2];
    link._odata[ss]()()(i1, 0) = subgroup._odata[ss]()()(1, 0) * link_elem[0][0] + subgroup._odata[ss]()()(1, 1) * link_elem[1][0];
    link._odata[ss]()()(i1, 1) = subgroup._odata[ss]()()(1, 0) * link_elem[0][1] + subgroup._odata[ss]()()(1, 1) * link_elem[1][1];
    link._odata[ss]()()(i1, 2) = subgroup._odata[ss]()()(1, 0) * link_elem[0][2] + subgroup._odata[ss]()()(1, 1) * link_elem[1][2];
  }
}

class i_Sigmas {
public:
  SU3::SU2Matrix ident;
  SU3::SU2Matrix pauli1;
  SU3::SU2Matrix pauli2;
  SU3::SU2Matrix pauli3;
  i_Sigmas() {
    ident = Complex(1.0);
    SU<2>::generator(0, pauli1);
    SU<2>::generator(1, pauli2);
    SU<2>::generator(2, pauli3);
    pauli1 = timesI(pauli1) * 2.0;
    pauli2 = timesI(pauli2) * 2.0;
    pauli3 = timesI(pauli3) * 2.0;
  }
};

void GF_SubGroupHeatBath(
       GridSerialRNG &sRNG, GridParallelRNG &pRNG,
       RealD coeff,  // coeff multiplying Re Tr(field * staple) in action; for Wilson: beta / 3.0; for S_GF1: beta * M^2
       LatticeColourMatrix &link,
       const LatticeColourMatrix &staple,  // multiplied by action coeffs so th
       int su2_subgroup, int cb, const std::string &table_path) {
     GridBase *rbGrid = link._grid;

     // static Integral_table integral_table("/home/yz/GFFA/jupyter/numerical_integration/lookup_table_M4_beta0.7796");
     static Integral_table integral_table(coeff, table_path);

     // Subgroup manipulation in the lie algebra space
     SU3::LatticeSU2Matrix u(rbGrid);  u.checkerboard = cb;// Kennedy pendleton "u" real projected normalised Sigma
     SU3::LatticeSU2Matrix ua(rbGrid);  ua.checkerboard = cb;// a in pauli form
     LatticeComplex udet(rbGrid);  udet.checkerboard = cb;// determinant of real(staple)

     my_su2Extract(udet, u, link, staple, su2_subgroup);

     // // FIXME: maybe this is necessary
     // // from the book:In the rare case that det A vanishes, any random link variable is accepted
     // LatticeComplex cone(grid); cone.checkerboard = cb;
     // LatticeReal adet(grid); adet.checkerboard = cb;
     // SU3::LatticeSU2Matrix lident(grid); lident.checkerboard = cb;
     // adet = abs(toReal(udet));
     // lident = Complex(1.0);
     // cone = Complex(1.0);
     // Real machine_epsilon = 1.0e-7;
     //
     // LatticeInteger tmp_int(grid); tmp_int.checkerboard = cb;
     // tmp_int = adet > machine_epsilon;
     // u = where(tmp_int, u, lident);
     // udet = where(tmp_int, udet, cone);

     LatticeComplex sqrt_udet(rbGrid); sqrt_udet.checkerboard = cb;
     sqrt_udet = sqrt(udet);

     LatticeReal k(rbGrid); k.checkerboard = cb;
     k = toReal(sqrt_udet); // FIXME: Wilson only; k = \sqrt{\det[staple]}

     std::vector<LatticeReal> a(4, rbGrid); for(auto &x: a) x.checkerboard = cb;

     LatticeReal tmp(rbGrid); tmp.checkerboard = cb;
     random(pRNG, tmp);
     parallel_for(int ss=0; ss<tmp._grid->oSites(); ++ss){ // ! cannot use grid->lSites() because of simd; use oSites()
       double *tmp_ptr = (double *)&tmp[ss];
       double *a0_ptr = (double *)&a[0][ss];
       double *k_ptr = (double *)&k[ss];
       for(int idx=0; idx<vReal::Nsimd(); idx+=2) {
         *(a0_ptr + idx) = integral_table.get_a0(*(k_ptr + idx), *(tmp_ptr + idx));
         *(a0_ptr + idx + 1) = *(a0_ptr + idx);
       }
     }

     //////////////////////////////////////////
     //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
     //////////////////////////////////////////
     LatticeReal a123mag(rbGrid); a123mag.checkerboard = cb;
     a123mag = sqrt(abs(1.0 - a[0] * a[0]));

     LatticeReal cos_theta(rbGrid); cos_theta.checkerboard = cb;
     LatticeReal sin_theta(rbGrid); sin_theta.checkerboard = cb;
     LatticeReal phi(rbGrid); phi.checkerboard = cb;

     random(pRNG, phi);
     const RealD twopi = 2.0 * M_PI;
     phi = phi * twopi;  // uniform in [0,2pi]
     random(pRNG, cos_theta);
     cos_theta = (cos_theta * 2.0) - 1.0;  // uniform in [-1,1]
     sin_theta = sqrt(abs(1.0 - cos_theta * cos_theta));

     a[1] = a123mag * sin_theta * cos(phi);
     a[2] = a123mag * sin_theta * sin(phi);
     a[3] = a123mag * cos_theta;

     static i_Sigmas i_sigmas;

     ua = toComplex(a[0]) * i_sigmas.ident + toComplex(a[1]) * i_sigmas.pauli1 +
          toComplex(a[2]) * i_sigmas.pauli2 + toComplex(a[3]) * i_sigmas.pauli3;

     // u = u * pow(sqrt_udet, -1.0);  //  u = link * staple / sqrt(det(staple)); // u is old X
     SU3::LatticeSU2Matrix uinv(rbGrid); uinv.checkerboard = cb;
     uinv = adj(u * pow(sqrt_udet, -1.0));

     SU3::LatticeSU2Matrix new_su2(rbGrid);   new_su2.checkerboard = cb;// rotated matrix after hb
     new_su2 = uinv * ua; // new su2 can be both uinv * ua or ua * uinv; they are both the same distribution.

     my_su2Insert(new_su2, link, su2_subgroup);


   }



}
}
