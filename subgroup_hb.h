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
  GridBase *grid(link.Grid());

  static int subgroup_index[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  int i0 = subgroup_index[su2_index][0];
  int i1 = subgroup_index[su2_index][1];

  autoView(subgroup_v, subgroup, AcceleratorWrite);  // It is fine to also read from this View
  // autoView(subgroup_v_read, subgroup, AcceleratorRead);  
  autoView(Determinant_v, Determinant, AcceleratorWrite);
  autoView(link_v, link, AcceleratorRead);
  autoView(staple_v, staple, AcceleratorRead);

  // thread_for(ss, grid->oSites(), {
  accelerator_for(ss, grid->oSites(), vcplx::Nsimd(), {
    subgroup_v[ss]()()(0, 0) = matrix_mult_elem(link_v[ss], staple_v[ss], i0, i0);
    subgroup_v[ss]()()(0, 1) = matrix_mult_elem(link_v[ss], staple_v[ss], i0, i1);
    subgroup_v[ss]()()(1, 0) = matrix_mult_elem(link_v[ss], staple_v[ss], i1, i0);
    subgroup_v[ss]()()(1, 1) = matrix_mult_elem(link_v[ss], staple_v[ss], i1, i1);

    subgroup_v[ss] = 0.5 * (subgroup_v[ss] - adj(subgroup_v[ss]) + trace(adj(subgroup_v[ss])));

    Determinant_v[ss] = subgroup_v[ss]()()(0, 0) * subgroup_v[ss]()()(1, 1)
                            - subgroup_v[ss]()()(0, 1) * subgroup_v[ss]()()(1, 0);
  });
}

template <class vcplx>
void my_su2Insert(const Lattice<SU3::iSU2Matrix<vcplx> > &subgroup,
                      Lattice<SU3::iSUnMatrix<vcplx> > &link, int su2_index) {
  GridBase *grid(link.Grid());

  static int subgroup_index[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  int i0 = subgroup_index[su2_index][0];
  int i1 = subgroup_index[su2_index][1];

  autoView(link_v, link, AcceleratorWrite);
  autoView(subgroup_v, subgroup, AcceleratorRead);

  // thread_for(ss, grid->oSites(), {
  accelerator_for(ss, grid->oSites(), vcplx::Nsimd(), {
    vcplx link_elem[2][3];
    for(int i=0; i<3; ++i) {
      link_elem[0][i] = link_v[ss]()()(i0, i); 
      link_elem[1][i] = link_v[ss]()()(i1, i);
    }

    //FIXME: Coalesce read/write

    link_v[ss]()()(i0, 0) = subgroup_v[ss]()()(0, 0) * link_elem[0][0] + subgroup_v[ss]()()(0, 1) * link_elem[1][0];
    link_v[ss]()()(i0, 1) = subgroup_v[ss]()()(0, 0) * link_elem[0][1] + subgroup_v[ss]()()(0, 1) * link_elem[1][1];
    link_v[ss]()()(i0, 2) = subgroup_v[ss]()()(0, 0) * link_elem[0][2] + subgroup_v[ss]()()(0, 1) * link_elem[1][2];
    link_v[ss]()()(i1, 0) = subgroup_v[ss]()()(1, 0) * link_elem[0][0] + subgroup_v[ss]()()(1, 1) * link_elem[1][0];
    link_v[ss]()()(i1, 1) = subgroup_v[ss]()()(1, 0) * link_elem[0][1] + subgroup_v[ss]()()(1, 1) * link_elem[1][1];
    link_v[ss]()()(i1, 2) = subgroup_v[ss]()()(1, 0) * link_elem[0][2] + subgroup_v[ss]()()(1, 1) * link_elem[1][2];
  });
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

LatticeComplex invCplx(const LatticeComplex& in) {

	LatticeComplex ret(in.Grid());
	ret.Checkerboard() = in.Checkerboard();

  autoView(ret_v, ret, AcceleratorWrite);
  autoView(in_v, in, AcceleratorRead);

	typename std::remove_const<typename std::remove_reference<decltype(in_v[0])>::type>::type one;
	one = 1.0;
  // thread_for(ss, in.Grid()->oSites(), {
  accelerator_for(ss, in.Grid()->oSites(), vComplex::Nsimd(), {
		ret_v[ss] = one / in_v[ss];
	});
	return ret;
}

void GF_SubGroupHeatBath(
       GridSerialRNG &sRNG, GridParallelRNG &pRNG,
       RealD coeff,  // coeff multiplying Re Tr(field * staple) in action; for Wilson: beta / 3.0; for S_GF1: beta * M^2
       LatticeColourMatrix &link,
       const LatticeColourMatrix &staple,  // multiplied by action coeffs so th
       int su2_subgroup, int cb, const std::string &table_path) {
     GridBase *rbGrid = link.Grid();

     static Integral_table integral_table(coeff, table_path);

     // Subgroup manipulation in the lie algebra space
     SU3::LatticeSU2Matrix u(rbGrid);  u.Checkerboard() = cb;// Kennedy pendleton "u" real projected normalised Sigma
     SU3::LatticeSU2Matrix ua(rbGrid);  ua.Checkerboard() = cb;// a in pauli form
     LatticeComplex udet(rbGrid);  udet.Checkerboard() = cb;// determinant of real(staple)

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

     LatticeComplex sqrt_udet(rbGrid); sqrt_udet.Checkerboard() = cb;
     sqrt_udet = sqrt(udet);

     LatticeReal k(rbGrid); k.Checkerboard() = cb;
     k = toReal(sqrt_udet); // FIXME: Wilson only; k = \sqrt{\det[staple]}

     std::vector<LatticeReal> a(4, rbGrid); for(auto &x: a) x.Checkerboard() = cb;

     LatticeReal tmp(rbGrid); tmp.Checkerboard() = cb;
     random(pRNG, tmp);

     autoView(a0_v, a[0], AcceleratorWrite);  // FIXME: we are both writing to and reading from a0_v
     autoView(tmp_v, tmp, AcceleratorRead);
     autoView(k_v, k, AcceleratorRead);


     // thread_for(ss, tmp.Grid()->oSites(), {
     accelerator_for(ss, tmp.Grid()->oSites(), vComplex::Nsimd(), {
       double *tmp_ptr = (double *)&tmp_v[ss];
       double *a0_ptr = (double *)&a0_v[ss];
       double *k_ptr = (double *)&k_v[ss];
       for(int idx=0; idx<vReal::Nsimd(); idx+=2) {
         *(a0_ptr + idx) = integral_table.get_a0(*(k_ptr + idx), *(tmp_ptr + idx));
         *(a0_ptr + idx + 1) = *(a0_ptr + idx);
       }
     });

     //////////////////////////////////////////
     //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
     //////////////////////////////////////////
     LatticeReal a123mag(rbGrid); a123mag.Checkerboard() = cb;
     // a123mag = sqrt(abs(1.0 - a[0] * a[0]));
     a123mag = sqrt(1.0 - a[0] * a[0]);

     LatticeReal cos_theta(rbGrid); cos_theta.Checkerboard() = cb;
     LatticeReal sin_theta(rbGrid); sin_theta.Checkerboard() = cb;
     LatticeReal phi(rbGrid); phi.Checkerboard() = cb;

     random(pRNG, phi);
     const RealD twopi = 2.0 * M_PI;
     phi = phi * twopi;  // uniform in [0,2pi]
     random(pRNG, cos_theta);
     cos_theta = (cos_theta * 2.0) - 1.0;  // uniform in [-1,1]
     // sin_theta = sqrt(abs(1.0 - cos_theta * cos_theta));
     sin_theta = sqrt(1.0 - cos_theta * cos_theta);

     a[1] = a123mag * sin_theta * cos(phi);
     a[2] = a123mag * sin_theta * sin(phi);
     a[3] = a123mag * cos_theta;

     static i_Sigmas i_sigmas;

     ua = toComplex(a[0]) * i_sigmas.ident + toComplex(a[1]) * i_sigmas.pauli1 +
          toComplex(a[2]) * i_sigmas.pauli2 + toComplex(a[3]) * i_sigmas.pauli3;

     SU3::LatticeSU2Matrix uinv(rbGrid); uinv.Checkerboard() = cb;
     // uinv = adj(u * pow(sqrt_udet, -1.0));
     uinv = adj(u * invCplx(sqrt_udet) ); // pow(, -1.0) is slow to calcualte

     SU3::LatticeSU2Matrix new_su2(rbGrid);   new_su2.Checkerboard() = cb;// rotated matrix after hb
     new_su2 = uinv * ua; // new su2 can be both uinv * ua or ua * uinv; they are both the same distribution.

     my_su2Insert(new_su2, link, su2_subgroup);

   }



}
}
