namespace Grid{
namespace QCD{


LatticeColourMatrix multField(const LatticeGaugeField &R1, const LatticeGaugeField &R2)
{
  LatticeColourMatrix R1mu(R1._grid);
  LatticeColourMatrix R2mu(R1._grid);
  LatticeColourMatrix ret(R1._grid);
  ret = zero;
  for (int mu = 0; mu < Nd; mu++) {
    R1mu = peekLorentz(R1, mu);
    R2mu = peekLorentz(R2, mu);
    ret += R1mu * R2mu;
  }
  return ret;
}

Real Hp(const LatticeGaugeField &P, const Momenta_k &KK)
{
  FFT theFFT((Grid::GridCartesian *)P._grid);

  LatticeGaugeField Pk(P._grid);
  theFFT.FFT_all_dim(Pk, P, FFT::forward);
  Pk = Pk * (1.0 / std::sqrt(KK.vol));

  LatticeColourMatrix sinKExpDotPk(P._grid);
  sinKExpDotPk = KK.sinKPsExpDotP_func(Pk);

  Real ret=0;
  LatticeColourMatrix ret1(P._grid);
  LatticeColourMatrix ret2(P._grid);

  ret1 = multField(adj(Pk) * (KK.one/KK.sinKEpsilonSquare), Pk);
  ret2 = adj(sinKExpDotPk) * sinKExpDotPk * KK.Ck_D;
  ret = TensorRemove(sum(trace( ret1 + ret2))).real();

  return ret;
}



LatticeGaugeField dHdP(LatticeGaugeField &P, const Momenta_k &KK)
{
  LatticeGaugeField ret(P._grid);
  FFT theFFT((Grid::GridCartesian *)P._grid);

  LatticeGaugeField realP(P._grid);
  realP = -timesI(P);  //this is necessary, otherwise nan.

  LatticeGaugeField Pk(P._grid);
  theFFT.FFT_all_dim(Pk, realP, FFT::forward);
  // Pk = Pk * (1.0 / std::sqrt(KK.vol)); //not necessary

  LatticeColourMatrix sinKExpDotPk(P._grid);
  sinKExpDotPk = KK.sinKPsExpDotP_func(Pk);

  LatticeColourMatrix dHdP2mu(P._grid);
  for(int mu=0; mu<Nd; ++mu)
  {
    dHdP2mu = KK.sinKNgExp[mu] * KK.Ck_D * sinKExpDotPk;
    pokeLorentz(ret, dHdP2mu, mu);
  }

  ret = ret + KK.one / KK.sinKEpsilonSquare * Pk;
  theFFT.FFT_all_dim(ret, ret, FFT::backward);
  ret = timesI(ret); // because realP = -timesI(P);
  return ret;
}

}
}
