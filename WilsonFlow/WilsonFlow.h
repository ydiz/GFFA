#ifndef ZYD_WILSONFLOW_H
#define ZYD_WILSONFLOW_H

#include <Grid/Grid.h>

namespace Grid {
namespace QCD {

class WilsonFlow_para {
public:
  std::vector<int> lat;
  double step_size;
  double adaptiveErrorTolerance;
  int Nstep;

  int StartTrajectory;
  int EndTrajectory;
  int TrajectoryInterval;
  std::string inFilePrefix;

  bool doSmear;
  bool saveSmearField;
  std::string smearFieldFilePrefix;

  bool calculateTopoCharge;
  std::string topoChargeOutFile;
};


template <class Gimpl>
class MyWilsonFlow {

    RealD initial_epsilon, epsilon, taus, adaptiveErrorTolerance, maxTau; // the taus printed out is flow time after this step
    int Nstep;
    bool hasCompleted = false; // for adaptive

    WilsonGaugeAction<Gimpl> SG;

    void evolve_step_adaptive(typename Gimpl::GaugeField&, RealD&);
    void evolve_step_adaptive_fixed_tau(typename Gimpl::GaugeField&);
    void evolve_step(typename Gimpl::GaugeField&) const;

 public:
    INHERIT_GIMPL_TYPES(Gimpl)

    explicit MyWilsonFlow(RealD _epsilon, RealD _adaptiveErrorTolerance, double _maxTau, int _Nstep=0):
        initial_epsilon(_epsilon),
        epsilon(_epsilon),
        adaptiveErrorTolerance(_adaptiveErrorTolerance),
        maxTau(_maxTau),
        Nstep(_Nstep),
        hasCompleted(false),
        SG(WilsonGaugeAction<Gimpl>(3.0)) { // WilsonGaugeAction with beta 3.0
            assert(epsilon > 0.0);
    }

    void smear(GaugeField& out, const GaugeField& in) const;
    void smear_adaptive(GaugeField&, const GaugeField&);
    void smear_adaptive_fixed_tau(GaugeField&, const GaugeField&);
};

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
void MyWilsonFlow<Gimpl>::smear(GaugeField& out, const GaugeField& in) const {
    out = in;
    double tau, energy_density;
    std::cout << GridLogMessage << "[WilsonFlow] step: "
              << 0 << "; tau: " << 0 << "; E: " << energyDensity(in) << std::endl;
    for (unsigned int step = 1; step <= Nstep; step++) {
        evolve_step(out);
        tau = step * epsilon;
        energy_density = energyDensity(out);
        std::cout << GridLogMessage << "[WilsonFlow] step: "
                  << step << "; tau: " << tau << "; E: "
                  << energy_density << "; t^2 E: " << tau * tau * energy_density << std::endl;
        std::cout << timeSliceTopologicalCharge(out) << std::endl;
    }
}


template <class Gimpl>
void MyWilsonFlow<Gimpl>::evolve_step_adaptive(typename Gimpl::GaugeField &U, double &tSqauredE) {

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


    // double new_tSqauredE = energyDensityPlaquette(U); // t^2 * <E>
    double energy_density = energyDensity(U);
    double new_tSqauredE = taus * taus * energy_density;


    if(step > 200) {
      std::cout << "Failed to converge!!!!!!" << std::endl;

      std::cout << GridLogMessage << "[WilsonFlow] step: "
        << step << "; tau: " << taus << "; E: "
        << energy_density << "; t^2 E: " << new_tSqauredE << std::endl;

      hasCompleted = true;
      return; // if hasCompleted, only update U
    }

    if(hasCompleted) {
      std::cout << GridLogMessage << "[WilsonFlow] step: "
        << step << "; tau: " << taus << "; E: "
        << energy_density << "; t^2 E: " << new_tSqauredE << std::endl;
      return; // if hasCompleted, only update U
    }
    if(new_tSqauredE > 0.3) { // if tSqauredE > 0.3, go back and use linear interpolation
      U = U0;
      taus = taus - epsilon;
      epsilon = epsilon * (0.3 - tSqauredE) / (new_tSqauredE - tSqauredE);
      taus = taus + epsilon;
      hasCompleted = true;
      evolve_step_adaptive(U, tSqauredE);
      step = 0; // reset step to 0 for the next smear_adaptive();
      return;
    }

    tSqauredE = new_tSqauredE;
    // calculate new step size
    Gimpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
    GaugeField diffU = U - Uprime;
    RealD diff = 1.0 / 9.0 * std::sqrt(maxNorm(diffU)); // d = 1/N^2 max_{x,mu} \sqrt( || U - Uprime || )


    // if d > δ the integration step is repeated; taus is unchanged.
    double new_epsilon;
    new_epsilon = epsilon * 0.95 * std::pow(adaptiveErrorTolerance/diff,1./3.);
    if(diff < adaptiveErrorTolerance) taus += new_epsilon; // taus is flow time after next step
    else {
      taus -= epsilon;
      U = U0;
      taus += new_epsilon;
    }

    // adjust integration step
    // epsilon = epsilon * 0.95 * std::pow(adaptiveErrorTolerance/diff,1./3.);
    epsilon = new_epsilon;
}






template <class Gimpl>
void MyWilsonFlow<Gimpl>::evolve_step_adaptive_fixed_tau(typename Gimpl::GaugeField &U) {

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


    taus += epsilon;
    double energy_density = energyDensity(U);
    double new_tSqauredE = taus * taus * energy_density;

      std::cout << "tau: " << taus << "; E: "
        << energy_density << "; t^2 E: " << new_tSqauredE << std::endl;


    // calculate new step size
    Gimpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
    GaugeField diffU = U - Uprime;
    RealD diff = 1.0 / 9.0 * std::sqrt(maxNorm(diffU)); // d = 1/N^2 max_{x,mu} \sqrt( || U - Uprime || )

    // if d > δ the integration step is repeated; taus is unchanged.
    if(diff > adaptiveErrorTolerance) {
      taus -= epsilon;
      U = U0;
    }


    double new_epsilon;
    new_epsilon = epsilon * 0.95 * std::pow(adaptiveErrorTolerance/diff,1./3.);

    epsilon = new_epsilon;
}



template <class Gimpl>
void MyWilsonFlow<Gimpl>::smear_adaptive(GaugeField& out, const GaugeField& in){
  out = in;
  epsilon = initial_epsilon;

  std::cout << GridLogDebug << "[WilsonFlow] step: "  << 0
    << "; tau: " << 0 << "; E: " << energyDensity(in) << std::endl;

  taus = epsilon; // initial taus is flow time after first step
  // taus = 0;
  hasCompleted = false;
  double tSqauredE = 0.;

  do{
    evolve_step_adaptive(out, tSqauredE);
    if(hasCompleted) break;
  } while(true);
}




template <class Gimpl>
void MyWilsonFlow<Gimpl>::smear_adaptive_fixed_tau(GaugeField& out, const GaugeField& in){
  out = in;
  epsilon = initial_epsilon;

  std::cout << GridLogDebug << "[WilsonFlow] step: "  << 0
    << "; tau: " << 0 << "; E: " << energyDensity(in) << std::endl;

  taus = 0;

  int step = 0;
  do{
    std::cout << GridLogMessage << "[WilsonFlow] step: " << step << "; ";
    step++;

    evolve_step_adaptive_fixed_tau(out);
    if(taus + epsilon > maxTau) {
      std::cout << GridLogMessage << "[WilsonFlow] step: " << step << "; ";
      epsilon = maxTau - taus;
      evolve_step_adaptive_fixed_tau(out);
      break;
    }
  } while(true);

}




}  // namespace QCD
} // namespace Grid

#endif // WILSONFLOW_H
