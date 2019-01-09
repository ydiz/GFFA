
namespace Grid{
namespace QCD{

// FIXME: change name: grid -> rbGrid

void GF_SubGroupHeatBath(
       GridSerialRNG &sRNG, GridParallelRNG &pRNG,
       RealD coeff,  // coeff multiplying Re Tr(field * staple) in action; for Wilson: beta / 3.0; for S_GF1: beta * M^2
       LatticeColourMatrix &link,
       const LatticeColourMatrix &staple,  // multiplied by action coeffs so th
       int su2_subgroup, int cb, const std::string &table_path) {
     GridBase *grid = link._grid;

     // static Integral_table integral_table("/home/yz/GFFA/jupyter/numerical_integration/lookup_table_M4_beta0.7796");
     static Integral_table integral_table(coeff, table_path);

     LatticeColourMatrix V(grid); V.checkerboard = cb;
     V = link * staple;

     // Subgroup manipulation in the lie algebra space
     SU3::LatticeSU2Matrix u(grid);  u.checkerboard = cb;// Kennedy pendleton "u" real projected normalised Sigma
     SU3::LatticeSU2Matrix uinv(grid); uinv.checkerboard = cb;
     SU3::LatticeSU2Matrix ua(grid);  ua.checkerboard = cb;// a in pauli form
     SU3::LatticeSU2Matrix b(grid);   b.checkerboard = cb;// rotated matrix after hb

     LatticeComplex udet(grid);  udet.checkerboard = cb;// determinant of real(staple)

     // Real part of Pauli decomposition
     // Note a subgroup can project to zero in cold start

     SU3::su2Extract(udet, u, V, su2_subgroup); // u = link * staple; udet = det(link * staple) = det(staple)
     u = u * 0.5; // FIXME: Grid defines su2Extract in a strange way!!!!
     udet = udet * 0.25; // FIXME: Grid defines su2Extract in a strange way!!!!

     // std::cout << u[0]<< "\t" << u[23]<< "\t"<< u[56] << std::endl; assert(0);

     // FIXME: maybe this is necessary
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

     LatticeComplex sqrt_udet(grid); sqrt_udet.checkerboard = cb;
     sqrt_udet = sqrt(udet);

     u = u * pow(sqrt_udet, -1.0);  //  u = link * staple / sqrt(det(staple)); // u is old X
     uinv = adj(u);

     LatticeReal k(grid); k.checkerboard = cb;
     k = toReal(sqrt_udet); // FIXME: Wilson only; k = \sqrt{\det[staple]}

     // uncomment the following to print out average k.
     // static double k_sum = 0;
     // static double k_num = 0;
     // {
       // std::cout << k[0]<< "\t" << k[23]<< "\t"<< k[56] << "\t" << k[34]<< "\t"
       //  << k[3]<< "\t" << k[35]<< "\t"<< k[58] << "\t" << k[32]<< "\t" << std::endl;
       // printRB(k);

       // std::cout << "max: " << maxLattice(k) << " min: " << minLattice(k) << std::endl;

        // double tmp_sum = sum(k)()()();
        // std::cout << "average: " << tmp_sum  / k._grid->gSites()<< std::endl;

        // k_sum += tmp_sum;
        // k_num += k._grid->gSites();
        // std::cout << "overall average: " << k_sum / k_num << std::endl;
        // assert(0);
     // }

     // LatticeReal a0_old(grid); a0_old.checkerboard = cb;
     // a0_old = toReal(0.5 * trace(u)); // old a0

     std::vector<LatticeReal> a(4, grid); for(auto &x: a) x.checkerboard = cb;

     LatticeReal tmp(grid); tmp.checkerboard = cb;
     random(pRNG, tmp);
     parallel_for(int ss=0; ss<tmp._odata.size(); ++ss){ // ! cannot use grid->lSites() because of simd; use oSites()
       double *tmp_ptr = (double *)&tmp[ss];
       double *a0_ptr = (double *)&a[0][ss];
       double *k_ptr = (double *)&k[ss];
       for(int idx=0; idx<vReal::Nsimd(); idx+=2) {
         *(a0_ptr + idx) = integral_table.get_a0(*(k_ptr + idx), *(tmp_ptr + idx));
         *(a0_ptr + idx + 1) = *(a0_ptr + idx);
       }
     }
      // std::cout << tmp[0]<< "\t" << a[0][0] <<  "\n"<< tmp[23]<< "\t"<< a[0][23] << "\n"<< tmp[46]<< "\t"<< a[0][46] <<std::endl;
      // assert(0);

     // LatticeReal metropolis_random(grid); metropolis_random.checkerboard = cb;
     // random(pRNG, metropolis_random);
     //
     // LatticeReal thresh(grid); thresh.checkerboard = cb;
     // thresh = exp(2.0 * coeff * (k - k_bar) * (a[0] - a0_old));
     //
     // LatticeInteger Accepted(grid); Accepted.checkerboard = cb;
     // Accepted = metropolis_random < thresh;
     //
     // // uncomment this to print out number of Accepted sites
     // {
     //   double numAccepted;
     //   LatticeReal rones(grid); rones.checkerboard = cb;
     //   rones = 1.0;
     //   LatticeReal rzeros(grid);rzeros.checkerboard = cb;
     //   rzeros = zero;
     //   LatticeReal rtmp(grid); rtmp.checkerboard = cb;
     //   rtmp = where(Accepted, rones, rzeros);
     //   numAccepted = sum(rtmp);
     //   std::cout << "numAccepted: " << numAccepted << std::endl;
     // }

     //////////////////////////////////////////
     //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
     //////////////////////////////////////////
     LatticeReal a123mag(grid); a123mag.checkerboard = cb;
     a123mag = sqrt(abs(1.0 - a[0] * a[0]));

     LatticeReal cos_theta(grid); cos_theta.checkerboard = cb;
     LatticeReal sin_theta(grid); sin_theta.checkerboard = cb;
     LatticeReal phi(grid); phi.checkerboard = cb;

     random(pRNG, phi);
     const RealD twopi = 2.0 * M_PI;
     phi = phi * twopi;  // uniform in [0,2pi]
     random(pRNG, cos_theta);
     cos_theta = (cos_theta * 2.0) - 1.0;  // uniform in [-1,1]
     sin_theta = sqrt(abs(1.0 - cos_theta * cos_theta));

     a[1] = a123mag * sin_theta * cos(phi);
     a[2] = a123mag * sin_theta * sin(phi);
     a[3] = a123mag * cos_theta;

     // these can be collected and moved to a static object.
     SU3::SU2Matrix ident = Complex(1.0);
     SU3::SU2Matrix pauli1;
     SU<2>::generator(0, pauli1);
     SU3::SU2Matrix pauli2;
     SU<2>::generator(1, pauli2);
     SU3::SU2Matrix pauli3;
     SU<2>::generator(2, pauli3);
     pauli1 = timesI(pauli1) * 2.0;
     pauli2 = timesI(pauli2) * 2.0;
     pauli3 = timesI(pauli3) * 2.0;
     // std::cout << pauli1 << std::endl;
     // std::cout << pauli2 << std::endl;
     // std::cout << pauli3 << std::endl;
     // assert(0);

     ua = toComplex(a[0]) * ident + toComplex(a[1]) * pauli1 +
          toComplex(a[2]) * pauli2 + toComplex(a[3]) * pauli3;

     b = uinv * ua;
     SU3::su2Insert(b, V, su2_subgroup);

     // link = where(Accepted, V * link, link);
     link = V * link;

   }



}
}
