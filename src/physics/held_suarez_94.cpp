#include <cmath>  // sin, cos, pow
#include "../parameter_input.hpp"
#include "../physics/held_suarez_94.hpp"
#include "../athena_math.hpp" // _sqr
#include "../globals.hpp" // Rgas

//! \file held_suarez_94.cpp
//  \brief implementation of held_suarez_94.hpp

HeldSuarez94::HeldSuarez94() {}

HeldSuarez94::HeldSuarez94(ParameterInput *pin)
{
  tdy = pin->GetReal("problem", "tdy");
  tdz = pin->GetReal("problem", "tdz");
  psrf = pin->GetReal("problem", "psrf");
  tsrf = pin->GetReal("problem", "tsrf");
  tmin = pin->GetReal("problem", "tmin");

  kappa = (pin->GetReal("hydro", "gamma") - 1.) / pin->GetReal("hydro", "gamma");
  rgas = Globals::Rgas / pin->GetReal("hydro", "mu");
  //rgas = Globals::my_rank / pin->GetReal("hydro", "mu");
  grav = - pin->GetReal("hydro", "grav_acc1");
}

Real HeldSuarez94::get_temp_eq(Real lat, Real pres)
{
  Real temp = (tsrf - tdy * _sqr(sin(lat)) - tdz * log(pres/psrf) * _sqr(cos(lat))) * pow(pres/psrf, kappa);
  return _max(tmin, temp);
}

Real HeldSuarez94::operator()(Real ptop)
{
  Real temp1 = get_temp_eq(lat, pbot),
       temp2 = get_temp_eq(lat, ptop);

  Real rho1 = pbot / (rgas * temp1),
       rho2 = ptop / (rgas * temp2);

  return (ptop - pbot)/dz + sqrt(rho1 * rho2) * grav;
}
