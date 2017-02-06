#ifndef HELD_SUAREZ_94_HPP
#define HELD_SUAREZ_94_HPP

#include "../athena.hpp"  // Real, MeshBlock

//! \file held_suarez_94.hpp
//  \brief classes and functions related to Held & Suarez, 94

// Forward declarations
class ParameterInput;

//! \class HeldSuarez94
//  \brief contains all Held & Suarez 94 parameters
class HeldSuarez94
{
  friend class MeshBlock;
public:
  HeldSuarez94();
  HeldSuarez94(ParameterInput *pin);
  Real get_temp_eq(Real lat, Real pres);
  Real operator()(Real ptop);

protected:
  Real tdy;
  Real tdz;
  Real tsrf;
  Real tmin;
  Real psrf;
  Real kappa;
  Real rgas;
  Real grav;

  Real lat;   //!< latitude
  Real pbot;  //!< pressure of the bottom cell
  Real dz;    //!< distance between two cells
};


#endif
