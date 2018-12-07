namespace Grid{
namespace QCD{

void GF_SubGroupHeatBath(
       GridSerialRNG &sRNG, GridParallelRNG &pRNG,
       RealD beta,  // coeff multiplying staple in action (with no 1/Nc)
       LatticeColourMatrix &link,
       const LatticeColourMatrix &barestaple,  // multiplied by action coeffs so th
       int su2_subgroup, int nheatbath, LatticeInteger &wheremask) {
     GridBase *grid = link._grid;

     int ntrials = 0;
     int nfails = 0;
     const RealD twopi = 2.0 * M_PI;

     LatticeColourMatrix staple(grid);

     staple = barestaple * (beta / 3.0);

     LatticeColourMatrix V(grid);
     V = link * staple;

     // Subgroup manipulation in the lie algebra space
     SU3::LatticeSU2Matrix u(
         grid);  // Kennedy pendleton "u" real projected normalised Sigma
     SU3::LatticeSU2Matrix uinv(grid);
     SU3::LatticeSU2Matrix ua(grid);  // a in pauli form
     SU3::LatticeSU2Matrix b(grid);   // rotated matrix after hb

     // Some handy constant fields
     LatticeComplex ones(grid);
     ones = 1.0;
     LatticeComplex zeros(grid);
     zeros = zero;
     LatticeReal rones(grid);
     rones = 1.0;
     LatticeReal rzeros(grid);
     rzeros = zero;
     LatticeComplex udet(grid);  // determinant of real(staple)
     LatticeInteger mask_true(grid);
     mask_true = 1;
     LatticeInteger mask_false(grid);
     mask_false = 0;

     /*
   PLB 156 P393 (1985) (Kennedy and Pendleton)
   Note: absorb "beta" into the def of sigma compared to KP paper; staple
   passed to this routine has "beta" already multiplied in
   Action linear in links h and of form:
       beta S = beta  Sum_p (1 - 1/Nc Re Tr Plaq )
   Writing Sigma = 1/Nc (beta Sigma') where sum over staples is "Sigma' "
        beta S = const - beta/Nc Re Tr h Sigma'
               = const - Re Tr h Sigma
   Decompose h and Sigma into (1, sigma_j) ; h_i real, h^2=1, Sigma_i complex
   arbitrary.
       Tr h Sigma = h_i Sigma_j Tr (sigma_i sigma_j)  = h_i Sigma_j 2 delta_ij
    Re Tr h Sigma = 2 h_j Re Sigma_j
   Normalised re Sigma_j = xi u_j
   With u_j a unit vector and U can be in SU(2);
   Re Tr h Sigma = 2 h_j Re Sigma_j = 2 xi (h.u)
   4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
    u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]
    xi = sqrt(Det)/2;
   Write a= u h in SU(2); a has pauli decomp a_j;
   Note: Product b' xi is unvariant because scaling Sigma leaves
         normalised vector "u" fixed; Can rescale Sigma so b' = 1.
     */

     ////////////////////////////////////////////////////////
     // Real part of Pauli decomposition
     // Note a subgroup can project to zero in cold start
     ////////////////////////////////////////////////////////
     SU3::su2Extract(udet, u, V, su2_subgroup);

     //////////////////////////////////////////////////////
     // Normalising this vector if possible; else identity
     //////////////////////////////////////////////////////
     LatticeComplex xi(grid);

     SU3::LatticeSU2Matrix lident(grid);

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

     LatticeComplex cone(grid);
     LatticeReal adet(grid);
     adet = abs(toReal(udet));
     lident = Complex(1.0);
     cone = Complex(1.0);
     Real machine_epsilon = 1.0e-7;
     u = where(adet > machine_epsilon, u, lident);
     udet = where(adet > machine_epsilon, udet, cone);

     xi = 0.5 * sqrt(udet);  // 4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
     u = 0.5 * u * pow(xi, -1.0);  //  u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

     // Debug test for sanity
     // uinv = adj(u);
     // b = u * uinv - 1.0;
     // assert(norm2(b) < 1.0e-4);

     /*
   Measure: Haar measure dh has d^4a delta(1-|a^2|)
   In polars:
     da = da0 r^2 sin theta dr dtheta dphi delta( 1 - r^2 -a0^2)
        = da0 r^2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r)(sqrt(1-a0^) +
   r) )
        = da0 r/2 sin theta dr dtheta dphi delta( (sqrt(1-a0^) - r) )
   Action factor Q(h) dh  = e^-S[h]  dh =  e^{  xi Tr uh} dh    // beta enters
   through xi
                                        =  e^{2 xi (h.u)} dh
                                        =  e^{2 xi h0u0}.e^{2 xi h1u1}.e^{2 xi
   h2u2}.e^{2 xi h3u3} dh
   Therefore for each site, take xi for that site
   i) generate  |a0|<1 with dist
      (1-a0^2)^0.5 e^{2 xi a0 } da0
   Take alpha = 2 xi  = 2 xi [ recall 2 beta/Nc unmod staple norm]; hence 2.0/Nc
   factor in Chroma ]
   A. Generate two uniformly distributed pseudo-random numbers R and R', R'',
   R''' in the unit interval;
   B. Set X = -(ln R)/alpha, X' =-(ln R')/alpha;
   C. Set C = cos^2(2pi R"), with R" another uniform random number in [0,1] ;
   D. Set A = XC;
   E. Let d  = X'+A;
   F. If R'''^2 :> 1 - 0.5 d,  go back to A;
   G. Set a0 = 1 - d;
   Note that in step D setting B ~ X - A and using B in place of A in step E will
   generate a second independent a 0 value.
     */

     /////////////////////////////////////////////////////////
     // count the number of sites by picking "1"'s out of hat
     /////////////////////////////////////////////////////////
     Integer hit = 0;
     LatticeReal rtmp(grid);
     rtmp = where(wheremask, rones, rzeros);
     RealD numSites = sum(rtmp);
     RealD numAccepted;
     LatticeInteger Accepted(grid);
     Accepted = zero;
     LatticeInteger newlyAccepted(grid);

     std::vector<LatticeReal> xr(4, grid);
     std::vector<LatticeReal> a(4, grid);
     LatticeReal d(grid);
     d = zero;
     LatticeReal alpha(grid);

     //    std::cout<<GridLogMessage<<"xi "<<xi <<std::endl;
     alpha = toReal(2.0 * xi);

     do {
       // A. Generate two uniformly distributed pseudo-random numbers R and R',
       // R'', R''' in the unit interval;
       random(pRNG, xr[0]);
       random(pRNG, xr[1]);
       random(pRNG, xr[2]);
       random(pRNG, xr[3]);

       // B. Set X = - ln R/alpha, X' = -ln R'/alpha
       xr[1] = -log(xr[1]) / alpha;
       xr[2] = -log(xr[2]) / alpha;

       // C. Set C = cos^2(2piR'')
       xr[3] = cos(xr[3] * twopi);
       xr[3] = xr[3] * xr[3];

       LatticeReal xrsq(grid);

       // D. Set A = XC;
       // E. Let d  = X'+A;
       xrsq = xr[2] + xr[1] * xr[3];

       d = where(Accepted, d, xr[2] + xr[1] * xr[3]);

       // F. If R'''^2 :> 1 - 0.5 d,  go back to A;
       LatticeReal thresh(grid);
       thresh = 1.0 - d * 0.5;
       xrsq = xr[0] * xr[0];
       LatticeInteger ione(grid);
       ione = 1;
       LatticeInteger izero(grid);
       izero = zero;

       newlyAccepted = where(xrsq < thresh, ione, izero);
       Accepted = where(newlyAccepted, newlyAccepted, Accepted);
       Accepted = where(wheremask, Accepted, izero);

       // FIXME need an iSum for integer to avoid overload on return type??
       rtmp = where(Accepted, rones, rzeros);
       numAccepted = sum(rtmp);

       hit++;

     } while ((numAccepted < numSites) && (hit < nheatbath));

     // G. Set a0 = 1 - d;
     a[0] = zero;
     a[0] = where(wheremask, 1.0 - d, a[0]);

     //////////////////////////////////////////
     //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
     //////////////////////////////////////////

     LatticeReal a123mag(grid);
     a123mag = sqrt(abs(1.0 - a[0] * a[0]));

     LatticeReal cos_theta(grid);
     LatticeReal sin_theta(grid);
     LatticeReal phi(grid);

     random(pRNG, phi);
     phi = phi * twopi;  // uniform in [0,2pi]
     random(pRNG, cos_theta);
     cos_theta = (cos_theta * 2.0) - 1.0;  // uniform in [-1,1]
     sin_theta = sqrt(abs(1.0 - cos_theta * cos_theta));

     a[1] = a123mag * sin_theta * cos(phi);
     a[2] = a123mag * sin_theta * sin(phi);
     a[3] = a123mag * cos_theta;

     ua = toComplex(a[0]) * ident + toComplex(a[1]) * pauli1 +
          toComplex(a[2]) * pauli2 + toComplex(a[3]) * pauli3;

     b = 1.0;
     b = where(wheremask, uinv * ua, b);
     SU3::su2Insert(b, V, su2_subgroup);

     // mask the assignment back based on Accptance
     link = where(Accepted, V * link, link);

     //////////////////////////////
     // Debug Checks
     // SU2 check
     // LatticeSU2Matrix check(grid);  // rotated matrix after hb
     // u = zero;
     // check = ua * adj(ua) - 1.0;
     // check = where(Accepted, check, u);
     // assert(norm2(check) < 1.0e-4);
     //
     // check = b * adj(b) - 1.0;
     // check = where(Accepted, check, u);
     // assert(norm2(check) < 1.0e-4);
     //
     // LatticeMatrix Vcheck(grid);
     // Vcheck = zero;
     // Vcheck = where(Accepted, V * adj(V) - 1.0, Vcheck);
     // //    std::cout<<GridLogMessage << "SU3 check " <<norm2(Vcheck)<<std::endl;
     // assert(norm2(Vcheck) < 1.0e-4);
     //
     // // Verify the link stays in SU(3)
     // //    std::cout<<GridLogMessage <<"Checking the modified link"<<std::endl;
     // Vcheck = link * adj(link) - 1.0;
     // assert(norm2(Vcheck) < 1.0e-4);
   /////////////////////////////////
   }



}
}
