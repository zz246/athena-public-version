//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hs94.cpp
//  \brief Problem generator for Held-Suarez-94 GCM bench mark.
//
// REFERENCE: I.M Held & M.J Suarez, "A Proposal for the Intercomparison of the
// Dynamical Cores of Atmospheric General Circulation Models"

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../athena_math.hpp" // _root
#include "../physics/held_suarez_94.hpp"  // HeldSuarez94

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // setup initial pressure/temperature field
  HeldSuarez94 hs(pin);
  double p1, t1;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      hs.lat = M_PI/2. - pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        phydro->w(IVX, i, j, k) = 0.;
        phydro->w(IVY, i, j, k) = 0.;
        phydro->w(IVZ, i, j, k) = 0.;
        if (i == is) {
          p1 = hs.psrf;
        } else {
          hs.dz = pcoord->x1v(i) - pcoord->x1v(i - 1);
          _root(hs.pbot, 0.5 * hs.pbot, 1., &p1, hs);
        }
        double t1 = hs.get_temp_eq(hs.lat, p1);
        phydro->w(IDN, i, j, k) = p1 / (hs.rgas * t1);
        phydro->w(IPR, i, j, k) = p1;
        hs.pbot = p1;
      }
    }

  // transfer to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
