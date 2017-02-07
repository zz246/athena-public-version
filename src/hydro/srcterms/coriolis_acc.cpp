//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to constant coriolis acceleration

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::Coriolis123
//  \brief Adds source terms for constant coriolis acceleration in axial direction to conserved variables

void HydroSourceTerms::Coriolis123(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  if (omega1_ != 0.0 || omega2_ != 0.0 || omega3_ != 0.0) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        cons(IM1,k,j,i) += 2.*dt*(omega3_*cons(IM2,k,j,i) - omega2_*cons(IM3,k,j,i));
        cons(IM2,k,j,i) += 2.*dt*(omega1_*cons(IM3,k,j,i) - omega3_*cons(IM1,k,j,i));
        cons(IM3,k,j,i) += 2.*dt*(omega2_*cons(IM1,k,j,i) - omega1_*cons(IM2,k,j,i));
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::CoriolisXYZ
//  \brief Adds source terms for constant coriolis acceleration in cartesian coordinate to conserved variables
void HydroSourceTerms::CoriolisXYZ(const Real dt,const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  Real omega1, omega2, omega3, theta, phi;

  if (omegax_ != 0.0 || omegay_ != 0.0 || omegaz_ != 0.0) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
#pragma omp parallel for schedule(static)
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        theta = pmb->pcoord->x2v(j);
        phi = pmb->pcoord->x3v(k);

        if (COORDINATE_SYSTEM == "cartesian") {
          omega1 = omegax_;
          omega2 = omegay_;
          omega3 = omegaz_;
        } else if (COORDINATE_SYSTEM == "cylindrical") {
          omega1 = cos(theta)*omegax_ + sin(theta)*omegay_;
          omega2 = - sin(theta)*omegax_ + cos(theta)*omegay_;
          omega3 = omegaz_;
        } else if (COORDINATE_SYSTEM == "spherical_polar") {
          omega1 = sin(theta)*cos(phi)*omegax_ + sin(theta)*sin(phi)*omegay_ + cos(theta)*omegaz_;
          omega2 = cos(theta)*cos(phi)*omegax_ + cos(theta)*sin(phi)*omegay_ - sin(theta)*omegaz_;
          omega3 = - sin(phi)*omegax_ + cos(phi)*omegay_;
        }

        cons(IM1,k,j,i) += 2.*dt*(omega3*cons(IM2,k,j,i) - omega2*cons(IM3,k,j,i));
        cons(IM2,k,j,i) += 2.*dt*(omega1*cons(IM3,k,j,i) - omega3*cons(IM1,k,j,i));
        cons(IM3,k,j,i) += 2.*dt*(omega2*cons(IM1,k,j,i) - omega1*cons(IM2,k,j,i));
      }
    }}
  }

  return;
}
