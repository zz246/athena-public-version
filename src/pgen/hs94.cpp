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

// functions for boundary conditions

void ProjectPressureInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// made global to share with BC function
static Real grav_acc;

//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // 3D problem
  // Enroll special BCs
  EnrollUserBoundaryFunction(INNER_X1, ProjectPressureInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, ProjectPressureOuterX1);
}

//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Held-Suarez problem generator
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // setup global variable shared by boundary condition
  grav_acc = - phydro->psrc->GetG1();

  // setup initial pressure/temperature field
  HeldSuarez94 hs(pin);
  Real p1, t1;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      hs.lat = M_PI/2. - pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        phydro->w(IVX, k, j, i) = 0.;
        phydro->w(IVY, k, j, i) = 0.;
        phydro->w(IVZ, k, j, i) = 0.;
        if (i == is) {
          p1 = hs.psrf;
        } else {
          hs.dz = pcoord->x1v(i) - pcoord->x1v(i - 1);
          int err = _root(hs.pbot, 0.5 * hs.pbot, 1., &p1, hs);
        }
        Real t1 = hs.get_temp_eq(hs.lat, p1);
        phydro->w(IDN, k, j, i) = p1 / (hs.rgas * t1);
        phydro->w(IPR, k, j, i) = p1;
        hs.pbot = p1;
      }
    }

  // transfer to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

//! \fn void ProjectPressureInnerX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
void ProjectPressureInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          if (n == IVX) {
            prim(IVX, k, j, is - i) = - prim(IVX, k, j, is + i - 1);  // reflect 1-vel
          } else if (n == IPR) {
            prim(IPR, k, j, is - i) = prim(IPR, k, j, is + i - 1)
              + prim(IDN, k, j, is + i - 1) * grav_acc * (2 * i - 1) * pco->dx1f(is + i - 1);
          } else {
            prim(n, k, j, is - i) = prim(n, k, j, is + i - 1);
          }
        }

  return;
}

//! \fn void ProjectPressureOuterX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = 1; i <= NGHOST; ++i) {
          if (n == IVX) {
            prim(IVX, k, j, ie + i) = - prim(IVX, k, j, ie - i + 1);  // reflect 1-vel
          } else if (n == IPR) {
            prim(IPR, k, j, ie + i) = prim(IPR, k, j, ie - i + 1)
              - prim(IDN, k, j, ie - i + 1) * grav_acc * (2 * i - 1) * pco->dx1f(ie - i + 1);
          } else {
            prim(n, k, j, ie + i) = prim(n, k, j, ie - i + 1);
          }
        }

  return;
}
