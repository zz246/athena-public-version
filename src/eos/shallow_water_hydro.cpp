//! \file shallow_water_hydro.cpp
//  \brief implements functions in class EquationOfState for shallow water hydrodynamics

// C/C++ headers
#include <cfloat>   // FLT_MIN

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"

// EquationOfState constructor
EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  density_floor_ = pin->GetOrAddReal("hydro", "dfloor", 1024 * FLT_MIN);
}

// destructor
EquationOfState::~EquationOfState()
{
}

void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
    AthenaArray<Real> const& prim_old, FaceField const& b,
    AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real di = 1. / cons(IDN,k,j,i);
        prim(IDN,k,j,i) = cons(IDN,k,j,i);
        prim(IVX,k,j,i) = cons(IM1,k,j,i) * di;
        prim(IVY,k,j,i) = cons(IM2,k,j,i) * di;
        prim(IVZ,k,j,i) = 0.;
      }
}

void EquationOfState::PrimitiveToConserved(AthenaArray<Real> const& prim,
    AthenaArray<Real> const& bc, AthenaArray<Real>& cons,
    Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real d = prim(IDN,k,j,i);
        cons(IDN,k,j,i) = d;
        cons(IM1,k,j,i) = prim(IVX,k,j,i) * d;
        cons(IM2,k,j,i) = prim(IVY,k,j,i) * d;
        cons(IM3,k,j,i) = 0.;
      }
}

Real EquationOfState::SoundSpeed(Real const prim[NHYDRO])
{
  return sqrt(prim[IDN]);
}
