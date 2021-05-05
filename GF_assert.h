namespace Grid{

bool isHermitian(const LatticeGaugeField &P)
{
  ColourMatrix M;
  peekSite(M, peekLorentz(P,0), std::vector<int>{1,2,3,4});

  return M()()(1,0) == conjugate(M()()(0,1)) && (M()()(0,0).imag() < 1e-7);
}
bool isHermitian(const LatticeColourMatrix &P)
{
  ColourMatrix M;
  peekSite(M, P, std::vector<int>{1,2,3,4});

  return M()()(1,0) == conjugate(M()()(0,1)) && (M()()(0,0).imag() < 1e-7);
}
bool isAntiHermitian(const LatticeGaugeField &P)
{
  ColourMatrix M;
  peekSite(M, peekLorentz(P,0), std::vector<int>{1,2,3,4});

  return M()()(1,0) == -conjugate(M()()(0,1)) && M()()(0,0) == -conjugate(M()()(0,0));
}
bool isAntiHermitian(const LatticeColourMatrix &P)
{
  ColourMatrix M;
  peekSite(M, P, std::vector<int>{1,2,3,4});

  return M()()(1,0) == -conjugate(M()()(0,1)) && M()()(0,0) == -conjugate(M()()(0,0));
}

ColourMatrix take00(const LatticeGaugeField &P)
{
  ColourMatrix m;
  peekSite(m, peekLorentz(P,0), std::vector<int>{0,0,0,0});
  return m;
}



}

