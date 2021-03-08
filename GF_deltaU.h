namespace Grid{
namespace QCD{


LatticeColourMatrix multField(const LatticeGaugeField &R1, const LatticeGaugeField &R2)
{
  LatticeColourMatrix R1mu(R1.Grid());
  LatticeColourMatrix R2mu(R1.Grid());
  LatticeColourMatrix ret(R1.Grid());
  ret = Zero();
  for (int mu = 0; mu < Nd; mu++) {
    R1mu = peekLorentz(R1, mu);
    R2mu = peekLorentz(R2, mu);
    ret += R1mu * R2mu;
  }
  return ret;
}

Real Hp(const LatticeGaugeField &P, const Momenta_k &KK)
{
  FFT theFFT((Grid::GridCartesian *)P.Grid());

  LatticeGaugeField Pk(P.Grid());
  theFFT.FFT_all_dim(Pk, P, FFT::forward);
  Pk = Pk * (1.0 / std::sqrt(KK.vol));

  LatticeColourMatrix sinKNgExpDotPk(P.Grid());
  sinKNgExpDotPk = KK.sinKNgExpDotP_func(Pk);

  Real ret=0;
  LatticeColourMatrix ret1(P.Grid());
  LatticeColourMatrix ret2(P.Grid());

  ret1 = multField(adj(Pk) * (KK.one/KK.FourSinKSquareEpsilon), Pk);
  ret2 = adj(sinKNgExpDotPk) * sinKNgExpDotPk * KK.Ck_D;
  ret = TensorRemove(sum(trace( ret1 + ret2))).real();

  return ret;
}



LatticeGaugeField dHdP(LatticeGaugeField &P, const Momenta_k &KK)
{
  // return P;  // FIXME
  GridBase *full_grid = P.Grid();
  GridBase *cell_grid = KK.Grid();
  Coordinate cell_size  = cell_grid->_fdimensions;

  LatticeGaugeField P_cell(cell_grid);
  localCopyRegion(P, P_cell, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);


  LatticeGaugeField ret_cell(cell_grid);
  FFT theFFT((Grid::GridCartesian *)cell_grid);

  LatticeGaugeField realP_cell(cell_grid);
  realP_cell = -timesI(P_cell);  //this is necessary, otherwise nan.

  LatticeGaugeField Pk_cell(cell_grid);
  theFFT.FFT_all_dim(Pk_cell, realP_cell, FFT::forward);
  // Pk = Pk * (1.0 / std::sqrt(KK.vol)); //not necessary

  LatticeColourMatrix sinKNgExpDotPk = KK.sinKNgExpDotP_func(Pk_cell);

  LatticeColourMatrix dHdP2mu(cell_grid);
  for(int mu=0; mu<Nd; ++mu)
  {
    dHdP2mu = KK.sinKPsExp[mu] * KK.Ck_D * sinKNgExpDotPk;
    pokeLorentz(ret_cell, dHdP2mu, mu);
  }

  ret_cell = ret_cell + KK.one / KK.FourSinKSquareEpsilon * Pk_cell;
  
  // Must do this when epsilon is 0 ;
  // set force of 0 mode to be 0 (KK.FourSinKSquareEpsilon = 0 for k=0)
  if(KK.epsi == 0.) {
    // std::cout << "I am setting the dH/dP of zero mode to 0; must do this when setting epsilon  to 0" << std::endl;
    typename LatticeGaugeField::vector_object::scalar_object m;
    m = 0.0;
    pokeSite(m, ret_cell, Coordinate({0,0,0,0}));
  }


  // std::cout << "Check dH/dP : " << std::endl;
  // // printGrid()_field_site(ret, {0,0,0,0});
  // print_grid_field_site(ret, {1,0,0,0});
  // print_grid_field_site(ret, {1,2,0,0});

  theFFT.FFT_all_dim(ret_cell, ret_cell, FFT::backward);
  ret_cell = timesI(ret_cell); // because realP before= -timesI(P);


  LatticeGaugeField ret(P.Grid()); ret = Zero();
  // LatticeGaugeField ret(P.Grid()); ret = P; // FIXME
  localCopyRegion(ret_cell, ret, Coordinate({0,0,0,0}), Coordinate({0,0,0,0}), cell_size);

  return ret;
}








}
}
