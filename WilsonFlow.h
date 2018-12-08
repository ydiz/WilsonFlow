#ifndef ZYD_WILSONFLOW_H
#define ZYD_WILSONFLOW_H

#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

template <class Gimpl>
class MyWilsonFlow: public Smear<Gimpl>{
    unsigned int Nstep;
    unsigned int measure_interval;
    mutable RealD epsilon, taus;

    mutable WilsonGaugeAction<Gimpl> SG;

    void evolve_step(typename Gimpl::GaugeField&) const;
    void evolve_step_adaptive(typename Gimpl::GaugeField&, RealD);
    RealD tau(unsigned int t)const {return epsilon*(t+1.0); }

 public:
    INHERIT_GIMPL_TYPES(Gimpl)

    explicit MyWilsonFlow(unsigned int Nstep, RealD epsilon, unsigned int interval = 1):
        Nstep(Nstep),
        epsilon(epsilon),
        measure_interval(interval),
        SG(WilsonGaugeAction<Gimpl>(3.0)) {
            // WilsonGaugeAction with beta 3.0
            assert(epsilon > 0.0);
            LogMessage();
    }

    void LogMessage() {
        std::cout << GridLogMessage
            << "[WilsonFlow] Nstep   : " << Nstep << std::endl;
        std::cout << GridLogMessage
            << "[WilsonFlow] epsilon : " << epsilon << std::endl;
        std::cout << GridLogMessage
            << "[WilsonFlow] full trajectory : " << Nstep * epsilon << std::endl;
    }

    virtual void smear(GaugeField&, const GaugeField&) const;

    virtual void derivative(GaugeField&, const GaugeField&, const GaugeField&) const {
        assert(0);
        // undefined for WilsonFlow
    }

    void smear_adaptive(GaugeField&, const GaugeField&, RealD maxTau);
    RealD energyDensityPlaquette(unsigned int step, const GaugeField& U) const;
    RealD energyDensityPlaquette(const GaugeField& U) const;
};


////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////
template <class Gimpl>
void MyWilsonFlow<Gimpl>::evolve_step(typename Gimpl::GaugeField &U) const{
    GaugeField Z(U._grid);
    GaugeField tmp(U._grid);
    SG.deriv(U, Z);
    Z *= 0.25;                                  // Z0 = 1/4 * F(U)
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0

    Z *= -17.0/8.0;
    SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
    Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
    Gimpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1

    Z *= -4.0/3.0;
    SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
    Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
    Gimpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
}

template <class Gimpl>
void MyWilsonFlow<Gimpl>::evolve_step_adaptive(typename Gimpl::GaugeField &U, RealD maxTau) {
    if (maxTau - taus < epsilon){
        epsilon = maxTau-taus;
    }
    //std::cout << GridLogMessage << "Integration epsilon : " << epsilon << std::endl;
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

    // Ramos
    Gimpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
    // Compute distance as norm^2 of the difference
    GaugeField diffU = U - Uprime;
    RealD diff = norm2(diffU);
    // adjust integration step
    //FIXME: zyd: should be diff = max_x,mu { sqrt(norm2(diffU(x, mu))) / N^2 }; do not know how to choose N
    // and  if d > δ the integration step is repeated.

    taus += epsilon;
    std::cout << GridLogMessage << "Adjusting integration step with distance: " << diff << std::endl;

    epsilon = epsilon*0.95*std::pow(1e-4/diff,1./3.);
    //std::cout << GridLogMessage << "New epsilon : " << epsilon << std::endl;

}

template <class Gimpl>
RealD MyWilsonFlow<Gimpl>::energyDensityPlaquette(unsigned int step, const GaugeField& U) const {
    RealD td = tau(step);
    return 2.0 * td * td * SG.S(U)/U._grid->gSites();
}

template <class Gimpl>
RealD MyWilsonFlow<Gimpl>::energyDensityPlaquette(const GaugeField& U) const {
    return 2.0 * taus * taus * SG.S(U)/U._grid->gSites();
}


template <class Gimpl>
void MyWilsonFlow<Gimpl>::smear(GaugeField& out, const GaugeField& in) const {
    out = in;
    for (unsigned int step = 1; step <= Nstep; step++) {
        auto start = std::chrono::high_resolution_clock::now();
        evolve_step(out);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        #ifdef WF_TIMING
        std::cout << "Time to evolve " << diff.count() << " s\n";
        #endif
        std::cout << GridLogMessage << "[WilsonFlow] Energy density (plaq) : "
		  << step << "  " << tau(step) << "  "
		  << energyDensityPlaquette(step,out) << std::endl;
         if( step % measure_interval == 0){
         std::cout << GridLogMessage << "[WilsonFlow] Top. charge           : "
            << step << "  "
            << WilsonLoops<PeriodicGimplR>::TopologicalCharge(out) << std::endl;
        }
    }
}

template <class Gimpl>
void MyWilsonFlow<Gimpl>::smear_adaptive(GaugeField& out, const GaugeField& in, RealD maxTau){
    out = in;
    taus = epsilon;
    unsigned int step = 0;
    do{
        step++;
        //std::cout << GridLogMessage << "Evolution time :"<< taus << std::endl;
        evolve_step_adaptive(out, maxTau);
        std::cout << GridLogMessage << "[WilsonFlow] Energy density (plaq) : "
		  << step << "  " << taus << "  "
		  << energyDensityPlaquette(out) << std::endl;
         if( step % measure_interval == 0){
         std::cout << GridLogMessage << "[WilsonFlow] Top. charge           : "
            << step << "  "
            << WilsonLoops<PeriodicGimplR>::TopologicalCharge(out) << std::endl;
        }
    } while (taus < maxTau);



}


}  // namespace QCD
} // namespace Grid

#endif // WILSONFLOW_H
