namespace Grid {

void stapleUpper(LatticeColourMatrix &staple, const LatticeGaugeField &Umu, int mu, int nu, int m=1, int n = 1) {
  assert(nu != mu);

  std::vector<LatticeColourMatrix> U(Nd, Umu.Grid());
  for (int d = 0; d < Nd; d++) {
    U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
  }

  LatticeColourMatrix tmp(Umu.Grid());
  tmp = PeriodicGimplR::CovShiftIdentityBackward(U[nu], nu);
  for(int i=0; i<n-1; ++i) tmp = PeriodicGimplR::CovShiftBackward(U[nu], nu, tmp);

  for(int i=0; i<m; ++i) tmp = PeriodicGimplR::CovShiftBackward(U[mu], mu, tmp);
  for(int i=0; i<n; ++i) tmp = PeriodicGimplR::CovShiftForward(U[nu], nu, tmp);

  staple = Cshift(tmp, mu, m);

  // staple = Cshift(
  //     PeriodicGimplR::CovShiftForward(
  //         U[nu], nu,
  //         PeriodicGimplR::CovShiftBackward(
  //             U[mu], mu, PeriodicGimplR::CovShiftIdentityBackward(U[nu], nu))),
  //     mu, 1);
  //
  // staple = Cshift(
  //     PeriodicGimplR::CovShiftForward(
  //         U[nu], nu, PeriodicGimplR::CovShiftForward(U[nu], nu,
  //         PeriodicGimplR::CovShiftBackward(
  //             U[mu], mu, PeriodicGimplR::CovShiftBackward(U[mu], mu, PeriodicGimplR::CovShiftBackward(U[nu], nu, PeriodicGimplR::CovShiftIdentityBackward(U[nu], nu))) ))),
  //     mu, 2);

}

////////////////////////////////////////////////////////////////////////
// the sum over all staples on each site in direction mu,nu, lower part
////////////////////////////////////////////////////////////////////////
void stapleLower(LatticeColourMatrix &staple, const LatticeGaugeField &Umu, int mu, int nu, int m=1, int n=1) {
  assert(nu != mu);

  std::vector<LatticeColourMatrix> U(Nd, Umu.Grid());
  for (int d = 0; d < Nd; d++) {
    U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
  }

  LatticeColourMatrix tmp(Umu.Grid());
  // tmp = PeriodicGimplR::CovShiftIdentityBackward(U[nu], nu));
  tmp = U[nu]; // p.s. PeriodicGimplR::CovShiftIdentityBackward does not change the matrix
  for(int i=0; i<n-1; ++i) tmp = PeriodicGimplR::CovShiftForward(U[nu], nu, tmp);

  for(int i=0; i<m; ++i) tmp = PeriodicGimplR::CovShiftBackward(U[mu], mu, tmp);
  for(int i=0; i<n; ++i) tmp = PeriodicGimplR::CovShiftBackward(U[nu], nu, tmp);

  staple = Cshift(tmp, mu, m);

  // staple = Cshift(
  //     PeriodicGimplR::CovShiftBackward(U[nu], nu,
  //                             PeriodicGimplR::CovShiftBackward(U[mu], mu, U[nu])),
  //     mu, 1);

  // staple = Cshift(
  //     PeriodicGimplR::CovShiftBackward(U[nu], nu, PeriodicGimplR::CovShiftBackward(U[nu], nu,
  //                             PeriodicGimplR::CovShiftBackward(U[mu], mu, PeriodicGimplR::CovShiftBackward(U[mu], mu,
  //                                   PeriodicGimplR::CovShiftForward(U[nu], nu, U[nu])) ))),
  //     mu, 2);
}


void fieldStrength(LatticeColourMatrix &FS, const LatticeGaugeField &Umu, int mu, int nu, int m=1, int n=1){

    LatticeColourMatrix Vup(Umu.Grid()), Vdn(Umu.Grid());
    stapleUpper(Vup, Umu, mu, nu, m, n);
    stapleLower(Vdn, Umu, mu, nu, m, n);
    LatticeColourMatrix v = Vup - Vdn;
    // LatticeColourMatrix u = PeekIndex<LorentzIndex>(Umu, mu);  // some redundant copies

    LatticeColourMatrix U_mu = PeekIndex<LorentzIndex>(Umu, mu);

    LatticeColourMatrix u = U_mu;
    for(int i=0; i<m-1; ++i) u = PeriodicGimplR::CovShiftForward(U_mu, mu, u);

    LatticeColourMatrix vu = v*u;
    //FS = 0.25*Ta(u*v + Cshift(vu, mu, -1));
    FS = (u*v + Cshift(vu, mu, -m)); //FS = (u*v + Cshift(vu, mu, -1));
    // LatticeColourMatrix tmp1 = adj(Vup * u);
    // LatticeColourMatrix tmp2 = Vdn * u;
    // FS = adj(u * Vup) - u * Vdn + Cshift(tmp1, mu, -m) - Cshift(tmp2, mu, -m);
    FS = 0.125*(FS - adj(FS));
}



std::vector<double> timeSliceTopologicalCharge_mn(const LatticeGaugeField &U, int m=1, int n=1) {
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  LatticeColourMatrix Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
  fieldStrength(Bx, U, Ydir, Zdir, m, n);
  fieldStrength(By, U, Zdir, Xdir, m, n);
  fieldStrength(Bz, U, Xdir, Ydir, m, n);

  LatticeColourMatrix tmp(U.Grid());
  if(m!=n) {
    fieldStrength(tmp, U, Ydir, Zdir, n, m);
    Bx = (Bx + tmp) * 0.5;
    fieldStrength(tmp, U, Zdir, Xdir, n, m);
    By = (By + tmp) * 0.5;
    fieldStrength(tmp, U, Xdir, Ydir, n, m);
    Bz = (Bz + tmp) * 0.5;
  }

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  LatticeColourMatrix Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
  fieldStrength(Ex, U, Tdir, Xdir, m, n);
  fieldStrength(Ey, U, Tdir, Ydir, m, n);
  fieldStrength(Ez, U, Tdir, Zdir, m, n);

  if(m!=n) {
    fieldStrength(tmp, U, Tdir, Xdir, n, m);
    Ex = (Ex + tmp) * 0.5;
    fieldStrength(tmp, U, Tdir, Ydir, n, m);
    Ey = (Ey + tmp) * 0.5;
    fieldStrength(tmp, U, Tdir, Zdir, n, m);
    Ez = (Ez + tmp) * 0.5;
  }

  double coeff = 8.0/(32.0*M_PI*M_PI) / m / m / n / n;

  LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);

  std::vector<typename LatticeComplex::vector_object::scalar_object> slice_sum;
  sliceSum(qfield, slice_sum, Tdir);
  std::vector<double> ret(slice_sum.size());
  for(int i=0; i<slice_sum.size(); ++i) ret[i] = TensorRemove(slice_sum[i]).real();

  // std::cout << m << ", " << n << ": " << ret << std::endl;
  return ret;
}

std::vector<double> timeSliceTopologicalCharge(const LatticeGaugeField &U) {
  double c1, c2, c3, c4, c5;

  c5 = 1.0 / 20.;
  c1 = (19. - 55. * c5) / 9.;
  c2 = (1. - 64. * c5) / 9.;
  c3 = (- 64. + 640. * c5) / 45.;
  c4 = 1./5. - 2 * c5;

  std::vector<double> ret(U.Grid()->_fdimensions[Tdir], 0);
  ret = c1 * timeSliceTopologicalCharge_mn(U, 1, 1) + c2 * timeSliceTopologicalCharge_mn(U, 2, 2)
        + c3 * timeSliceTopologicalCharge_mn(U, 1, 2) + c4 * timeSliceTopologicalCharge_mn(U, 1, 3)
        + c5 * timeSliceTopologicalCharge_mn(U, 3, 3);
  // ret += c3 * timeSliceTopologicalCharge_mn(U, 1, 2) + ;
  // ret += c4 * timeSliceTopologicalCharge_mn(U, 1, 3);
  // ret += c5 * timeSliceTopologicalCharge_mn(U, 3, 3);

  return ret;

}

std::vector<double> timeSliceTopologicalCharge_1Li(const LatticeGaugeField &U) {

  std::vector<double> ret(U.Grid()->_fdimensions[Tdir], 0);
  ret = timeSliceTopologicalCharge_mn(U, 1, 1);

  return ret;

}

double globalTopologicalCharge_1Li(const LatticeGaugeField &U){
  std::vector<double> timeSliceTC;
  timeSliceTC = timeSliceTopologicalCharge_1Li(U);
  double ret = 0;
  for(double x: timeSliceTC) ret += x;
  return ret;
}

double globalTopologicalCharge(const LatticeGaugeField &U){
  std::vector<double> timeSliceTC;
  timeSliceTC = timeSliceTopologicalCharge(U);
  double ret = 0;
  for(double x: timeSliceTC) ret += x;
  return ret;
}

// double globalTopologicalCharge(const LatticeGaugeField &U){
//     // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
//     LatticeColourMatrix Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
//     fieldStrength(Bx, U, Ydir, Zdir);
//     fieldStrength(By, U, Zdir, Xdir);
//     fieldStrength(Bz, U, Xdir, Ydir);
//
//     // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
//     LatticeColourMatrix Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
//     fieldStrength(Ex, U, Tdir, Xdir);
//     fieldStrength(Ey, U, Tdir, Ydir);
//     fieldStrength(Ez, U, Tdir, Zdir);
//
//     double coeff = 8.0/(32.0*M_PI*M_PI);
//
//     LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
//     auto Tq = sum(qfield);
//     return TensorRemove(Tq).real();
// }

double energyDensity_mn(const LatticeGaugeField &U, int m=1, int n=1) {
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  LatticeColourMatrix Bx(U.Grid()), By(U.Grid()), Bz(U.Grid());
  fieldStrength(Bx, U, Ydir, Zdir, m, n);
  fieldStrength(By, U, Zdir, Xdir, m, n);
  fieldStrength(Bz, U, Xdir, Ydir, m, n);

  LatticeColourMatrix tmp(U.Grid());
  if(m!=n) {
    fieldStrength(tmp, U, Ydir, Zdir, n, m);
    Bx = (Bx + tmp) * 0.5;
    fieldStrength(tmp, U, Zdir, Xdir, n, m);
    By = (By + tmp) * 0.5;
    fieldStrength(tmp, U, Xdir, Ydir, n, m);
    Bz = (Bz + tmp) * 0.5;
  }

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  LatticeColourMatrix Ex(U.Grid()), Ey(U.Grid()), Ez(U.Grid());
  fieldStrength(Ex, U, Tdir, Xdir, m, n);
  fieldStrength(Ey, U, Tdir, Ydir, m, n);
  fieldStrength(Ez, U, Tdir, Zdir, m, n);

  if(m!=n) {
    fieldStrength(tmp, U, Tdir, Xdir, n, m);
    Ex = (Ex + tmp) * 0.5;
    fieldStrength(tmp, U, Tdir, Ydir, n, m);
    Ey = (Ey + tmp) * 0.5;
    fieldStrength(tmp, U, Tdir, Zdir, n, m);
    Ez = (Ez + tmp) * 0.5;
  }

  // double coeff = 8.0/(32.0*M_PI*M_PI);
  // LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
  double coeff = 2.0 / 2.0;

  LatticeComplex qfield = coeff*trace(Bx*Bx + By*By + Bz*Bz + Ex*Ex + Ey*Ey + Ez*Ez);

  // std::vector<typename LatticeComplex::vector_object::scalar_object> slice_sum;
  // sliceSum(qfield, slice_sum, Tdir);
  // std::vector<double> ret(slice_sum.size());
  // for(int i=0; i<slice_sum.size(); ++i) ret[i] = TensorRemove(slice_sum[i]).real();
  std::complex<double> ttt = sum(qfield)()()() / double(U.Grid()->gSites());
  // std::cout << ttt << std::endl;
  double ret = ttt.real();
  // double ret = TensorRemove(tmp).real();
  return ret;
}

double energyDensity(const LatticeGaugeField &U) {
  // double c1, c2, c3, c4, c5;
  //
  // c5 = 1.0 / 20.;
  // c1 = (19. - 55. * c5) / 9.;
  // c2 = (1. - 64. * c5) / 9.;
  // c3 = (- 64. + 640. * c5) / 45.;
  // c4 = 1./5. - 2 * c5;
  //
  // std::vector<double> ret(U.Grid()->_fdimensions[Tdir], 0);
  // ret = c1 * energyDensity_mn(U, 1, 1) + c2 * energyDensity_mn(U, 2, 2)
  //       + c3 * energyDensity_mn(U, 1, 2) + c4 * energyDensity_mn(U, 1, 3)
  //       + c5 * energyDensity_mn(U, 3, 3);
  double ret;
  ret = - energyDensity_mn(U, 1, 1); //FIXME: do not know why multipying by -1 is necessary

  return ret;

}


}

