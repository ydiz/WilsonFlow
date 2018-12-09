#ifndef ZYD_WILSONFLOW_H
#define ZYD_WILSONFLOW_H

#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

double maxNorm(const LatticeGaugeField& U) {
  Lattice<iVector<iScalar<iScalar<vRealD>>, Nd> > U_norm2(U._grid);

  parallel_for(int ss=0;ss<U_norm2._grid->oSites();ss++){
    for(int mu=0; mu<Nd; mu++) {
      U_norm2[ss](mu)()() = 0.;
      for(int c1=0;c1<Nc;c1++)
        for(int c2=0;c2<Nc;c2++)
            U_norm2[ss](mu)()() += toReal(U[ss](mu)()(c1, c2) * U[ss](mu)()(c1, c2));
    }
  }

  typedef typename decltype(U_norm2)::scalar_type scalar_type;

  double max_val = 0.;
  #pragma omp parallel for reduction(max : max_val)
  for(int ss=0;ss<U_norm2._grid->oSites();ss++)
  {
    for(int mu=0; mu<Nd; mu++) {
      scalar_type *sobj = (scalar_type *)& U_norm2[ss](mu)()();
      for(int idx=0; idx<U_norm2._grid->Nsimd(); ++idx)
        if( *(sobj + idx) > max_val)
            max_val = *(sobj + idx);
    }
  }

  #ifndef GRID_COMMS_NONE
  MPI_Allreduce(MPI_IN_PLACE, &max_val, 1, MPI_DOUBLE, MPI_MAX, U_norm2._grid->communicator);
  #endif

  return max_val;
}


class WilsonFlow_para {
public:
  int steps;
  double step_size;
  int measure_interval;
  double maxTau;
  double adaptiveErrorTolerance;
};

template <class Gimpl>
class MyWilsonFlow {
    // unsigned int Nstep;
    // unsigned int measure_interval;
    mutable RealD epsilon, taus, adaptiveErrorTolerance; // the taus printed out is flow time after this step
    bool hasCompleted = false; // for adaptive

    mutable WilsonGaugeAction<Gimpl> SG;

    void evolve_step_adaptive(typename Gimpl::GaugeField&, RealD&);

 public:
    INHERIT_GIMPL_TYPES(Gimpl)

    explicit MyWilsonFlow(unsigned int Nstep, RealD epsilon, RealD _adaptiveErrorTolerance, unsigned int interval = 1):
        epsilon(epsilon),
        adaptiveErrorTolerance(_adaptiveErrorTolerance),
        hasCompleted(false),
        SG(WilsonGaugeAction<Gimpl>(3.0)) { // WilsonGaugeAction with beta 3.0
            assert(epsilon > 0.0);
    }

    void smear_adaptive(GaugeField&, const GaugeField&);
    RealD energyDensityPlaquette(const GaugeField& U) const;
};


template <class Gimpl>
void MyWilsonFlow<Gimpl>::evolve_step_adaptive(typename Gimpl::GaugeField &U, double &tSqauredAveE) {

    static int step = 0;
    step++;
    LatticeGaugeField U0 = U;
    // calculate Uflow (third order method) and Uflow^prime (second order method)
    GaugeField Z(U._grid);
    GaugeField Zprime(U._grid);
    GaugeField tmp(U._grid), Uprime(U._grid);
    Uprime = U;
    SG.deriv(U, Z);
    Zprime = -Z;
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Zprime += 2.0*tmp;                          // zyd: Zprime = - Z0 + 2 Z1
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    Gimpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2


    double new_tSqauredAveE = energyDensityPlaquette(U); // t^2 * <E>
    std::cout << GridLogMessage << "[WilsonFlow] Energy density (plaq) : "
              << step << "  " << taus << "  " << new_tSqauredAveE << std::endl;

    if(hasCompleted) return; // if hasCompleted, only update U
    if(new_tSqauredAveE > 0.3) { // if tSqauredAveE > 0.3, go back and use linear interpolation
      U = U0;
      taus = taus - epsilon;
      epsilon = epsilon * (0.3 - tSqauredAveE) / (new_tSqauredAveE - tSqauredAveE);
      taus = taus + epsilon;
      hasCompleted = true;
      evolve_step_adaptive(U, tSqauredAveE);
      step = 0; // reset step to 0 for the next smear_adaptive();
      return;
    }

    tSqauredAveE = new_tSqauredAveE;
    // calculate new step size
    Gimpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
    GaugeField diffU = U - Uprime;
    RealD diff = 1.0 / 9.0 * std::sqrt(maxNorm(diffU)); // d = 1/N^2 max_{x,mu} \sqrt( || U - Uprime || )
    // if d > δ the integration step is repeated; taus is unchanged.
    if(diff < adaptiveErrorTolerance) taus += epsilon;
    else {
      U = U0;
    }
    // std::cout << GridLogMessage << "Adjusting integration step with distance: " << diff << std::endl;
    // adjust integration step
    epsilon = epsilon * 0.95 * std::pow(adaptiveErrorTolerance/diff,1./3.);
}

template <class Gimpl>
RealD MyWilsonFlow<Gimpl>::energyDensityPlaquette(const GaugeField& U) const {
    return 2.0 * taus * taus * SG.S(U)/U._grid->gSites();
}


template <class Gimpl>
void MyWilsonFlow<Gimpl>::smear_adaptive(GaugeField& out, const GaugeField& in){
    out = in;
    taus = epsilon; // initial taus is flow time after first step
    // taus = 0; //zyd: I changed this to: initial taus = 0;
    // unsigned int step = 0;
    hasCompleted = false;
    double tSqauredAveE = 0.;
    do{
        // step++;
        //std::cout << GridLogMessage << "Evolution time :"<< taus << std::endl;
        evolve_step_adaptive(out, tSqauredAveE);

        if(hasCompleted) break;
        //   std::cout << GridLogMessage << "[WilsonFlow] Energy density (plaq) : "
  		  // << step << "  " << taus << "  " << energyDensityPlaquette(out) << std::endl;
    } while(true);
}


}  // namespace QCD
} // namespace Grid

#endif // WILSONFLOW_H
