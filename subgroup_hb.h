namespace Grid{
namespace QCD{

// FIXME: change name: grid -> rbGrid

void GF_SubGroupHeatBath(
       GridSerialRNG &sRNG, GridParallelRNG &pRNG,
       RealD beta,  // coeff multiplying staple in action (with no 1/Nc)
       LatticeColourMatrix &link,
       const LatticeColourMatrix &barestaple,  // multiplied by action coeffs so th
       int su2_subgroup, int nheatbath, int cb) {
     GridBase *grid = link._grid;

     int ntrials = 0;
     int nfails = 0;
     const RealD twopi = 2.0 * M_PI;

     LatticeColourMatrix staple(grid); staple.checkerboard = cb;

     staple = barestaple * (beta / 3.0);

     LatticeColourMatrix V(grid); V.checkerboard = cb;
     V = link * staple;

     // Subgroup manipulation in the lie algebra space
     SU3::LatticeSU2Matrix u(grid);  u.checkerboard = cb;// Kennedy pendleton "u" real projected normalised Sigma
     SU3::LatticeSU2Matrix uinv(grid); uinv.checkerboard = cb;
     SU3::LatticeSU2Matrix ua(grid);  ua.checkerboard = cb;// a in pauli form
     SU3::LatticeSU2Matrix b(grid);   b.checkerboard = cb;// rotated matrix after hb

     // Some handy constant fields
     LatticeReal rones(grid); rones.checkerboard = cb;
     rones = 1.0;
     LatticeReal rzeros(grid);rzeros.checkerboard = cb;
     rzeros = zero;
     LatticeComplex udet(grid);  udet.checkerboard = cb;// determinant of real(staple)

     ////////////////////////////////////////////////////////
     // Real part of Pauli decomposition
     // Note a subgroup can project to zero in cold start
     ////////////////////////////////////////////////////////
     SU3::su2Extract(udet, u, V, su2_subgroup);

     //////////////////////////////////////////////////////
     // Normalising this vector if possible; else identity
     //////////////////////////////////////////////////////
     LatticeComplex xi(grid); xi.checkerboard = cb;

     SU3::LatticeSU2Matrix lident(grid); lident.checkerboard = cb;

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

     LatticeComplex cone(grid); cone.checkerboard = cb;
     LatticeReal adet(grid); adet.checkerboard = cb;
     adet = abs(toReal(udet));
     lident = Complex(1.0);
     cone = Complex(1.0);
     Real machine_epsilon = 1.0e-7;

     LatticeInteger tmp_int(grid); tmp_int.checkerboard = cb;
     tmp_int = adet > machine_epsilon;
     u = where(tmp_int, u, lident);
     udet = where(tmp_int, udet, cone);


     xi = 0.5 * sqrt(udet);  // 4xi^2 = Det [ Sig - Sig^dag  + 1 Tr Sigdag]
     u = 0.5 * u * pow(xi, -1.0);  //  u   = 1/2xi [ Sig - Sig^dag  + 1 Tr Sigdag]

     uinv = adj(u);

     /////////////////////////////////////////////////////////
     // count the number of sites by picking "1"'s out of hat
     /////////////////////////////////////////////////////////
     Integer hit = 0;
     LatticeReal rtmp(grid); rtmp.checkerboard = cb;
     // rtmp = 1.0;

     RealD numSites = grid->gSites(); //sum(rtmp);
     RealD numAccepted;
     LatticeInteger Accepted(grid); Accepted.checkerboard = cb;
     Accepted = zero;
     LatticeInteger newlyAccepted(grid); newlyAccepted.checkerboard = cb;

     std::vector<LatticeReal> xr(4, grid); for(auto &x: xr) x.checkerboard = cb;
     std::vector<LatticeReal> a(4, grid); for(auto &x: a) x.checkerboard = cb;
     LatticeReal d(grid); d.checkerboard = cb;
     d = zero;
     LatticeReal alpha(grid); alpha.checkerboard = cb;

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
       LatticeReal thresh(grid); thresh.checkerboard = cb;
       thresh = 1.0 - d * 0.5;
       xrsq = xr[0] * xr[0];
       LatticeInteger ione(grid); ione.checkerboard = cb;
       ione = 1;
       LatticeInteger izero(grid); izero.checkerboard = cb;
       izero = zero;

       newlyAccepted = where(xrsq < thresh, ione, izero);
       Accepted = where(newlyAccepted, newlyAccepted, Accepted);

       // FIXME need an iSum for integer to avoid overload on return type??
       rtmp = where(Accepted, rones, rzeros);
       numAccepted = sum(rtmp);

       hit++;

     } while ((numAccepted < numSites) && (hit < nheatbath));

     // G. Set a0 = 1 - d;
     a[0] = 1.0 - d;

     //////////////////////////////////////////
     //    ii) generate a_i uniform on two sphere radius (1-a0^2)^0.5
     //////////////////////////////////////////

     LatticeReal a123mag(grid); a123mag.checkerboard = cb;
     a123mag = sqrt(abs(1.0 - a[0] * a[0]));

     LatticeReal cos_theta(grid); cos_theta.checkerboard = cb;
     LatticeReal sin_theta(grid); sin_theta.checkerboard = cb;
     LatticeReal phi(grid); phi.checkerboard = cb;

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

     b = uinv * ua;
     SU3::su2Insert(b, V, su2_subgroup);

     link = where(Accepted, V * link, link);

   }



}
}
