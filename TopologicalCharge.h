namespace Grid {
namespace QCD {

// fieldStrength -> fieldStrength

//FIXME: specify Gimpl

void stapleUpper(LatticeColourMatrix &staple, const LatticeGaugeField &Umu, int mu, int nu, int m=1, int n = 1) {
  assert(nu != mu);

  std::vector<LatticeColourMatrix> U(Nd, Umu._grid);
  for (int d = 0; d < Nd; d++) {
    U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
  }

  LatticeColourMatrix tmp(Umu._grid);
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

  std::vector<LatticeColourMatrix> U(Nd, Umu._grid);
  for (int d = 0; d < Nd; d++) {
    U[d] = PeekIndex<LorentzIndex>(Umu, d);// some redundant copies
  }

  LatticeColourMatrix tmp(Umu._grid);
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

    LatticeColourMatrix Vup(Umu._grid), Vdn(Umu._grid);
    stapleUpper(Vup, Umu, mu, nu, m, n);
    stapleLower(Vdn, Umu, mu, nu, m, n);
    LatticeColourMatrix v = Vup - Vdn;
    // LatticeColourMatrix u = PeekIndex<LorentzIndex>(Umu, mu);  // some redundant copies

    // LatticeColourMatrix u = ;
    LatticeColourMatrix U_mu = PeekIndex<LorentzIndex>(Umu, mu);

    LatticeColourMatrix u = U_mu;
    for(int i=0; i<m-1; ++i) u = PeriodicGimplR::CovShiftForward(U_mu, mu, u);

    LatticeColourMatrix vu = v*u;
    //FS = 0.25*Ta(u*v + Cshift(vu, mu, -1));
    FS = (u*v + Cshift(vu, mu, -m)); //FS = (u*v + Cshift(vu, mu, -1));
    FS = 0.125*(FS - adj(FS));
}



std::vector<double> timeSliceTopologicalCharge(const LatticeGaugeField &U, int m=1, int n=1) {
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  LatticeColourMatrix Bx(U._grid), By(U._grid), Bz(U._grid);
  fieldStrength(Bx, U, Ydir, Zdir, m, n);
  fieldStrength(By, U, Zdir, Xdir, m, n);
  fieldStrength(Bz, U, Xdir, Ydir, m, n);

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  LatticeColourMatrix Ex(U._grid), Ey(U._grid), Ez(U._grid);
  fieldStrength(Ex, U, Tdir, Xdir, m, n);
  fieldStrength(Ey, U, Tdir, Ydir, m, n);
  fieldStrength(Ez, U, Tdir, Zdir, m, n);

  double coeff = 8.0/(32.0*M_PI*M_PI);

  LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);

  std::vector<typename LatticeComplex::vector_object::scalar_object> slice_sum;
  sliceSum(qfield, slice_sum, Tdir);
  std::vector<double> ret(slice_sum.size());
  for(int i=0; i<slice_sum.size(); ++i) ret[i] = TensorRemove(slice_sum[i]).real();
  return ret;
}

// double globalTopologicalCharge(const LatticeGaugeField &U){
//     // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
//     LatticeColourMatrix Bx(U._grid), By(U._grid), Bz(U._grid);
//     fieldStrength(Bx, U, Ydir, Zdir);
//     fieldStrength(By, U, Zdir, Xdir);
//     fieldStrength(Bz, U, Xdir, Ydir);
//
//     // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
//     LatticeColourMatrix Ex(U._grid), Ey(U._grid), Ez(U._grid);
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

}
}
