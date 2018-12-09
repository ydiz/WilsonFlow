namespace Grid {
namespace QCD {

double globalTopologicalCharge(const LatticeGaugeField &U){
    // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
    LatticeColourMatrix Bx(U._grid), By(U._grid), Bz(U._grid);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, U, Ydir, Zdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(By, U, Zdir, Xdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, U, Xdir, Ydir);

    // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
    LatticeColourMatrix Ex(U._grid), Ey(U._grid), Ez(U._grid);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ex, U, Tdir, Xdir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ey, U, Tdir, Ydir);
    WilsonLoops<PeriodicGimplR>::FieldStrength(Ez, U, Tdir, Zdir);

    double coeff = 8.0/(32.0*M_PI*M_PI);

    LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);
    auto Tq = sum(qfield);
    return TensorRemove(Tq).real();
}

std::vector<double> timeSliceTopologicalCharge(const LatticeGaugeField &U) {
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  LatticeColourMatrix Bx(U._grid), By(U._grid), Bz(U._grid);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, U, Ydir, Zdir);
  WilsonLoops<PeriodicGimplR>::FieldStrength(By, U, Zdir, Xdir);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, U, Xdir, Ydir);

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  LatticeColourMatrix Ex(U._grid), Ey(U._grid), Ez(U._grid);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Ex, U, Tdir, Xdir);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Ey, U, Tdir, Ydir);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Ez, U, Tdir, Zdir);

  double coeff = 8.0/(32.0*M_PI*M_PI);

  LatticeComplex qfield = coeff*trace(Bx*Ex + By*Ey + Bz*Ez);

  std::vector<typename LatticeComplex::vector_object::scalar_object> slice_sum;
  sliceSum(qfield, slice_sum, Tdir);
  std::vector<double> ret(slice_sum.size());
  for(int i=0; i<slice_sum.size(); ++i) ret[i] = TensorRemove(slice_sum[i]).real();
  return ret;
}



}
}
