//! \file jupiter_pole.cpp
//  \brief jupiter polar model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../athena_math.hpp"     // _root, _sqr
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp" // ran2
#include "../particle/particle.hpp"

// support functions
void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0);
void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0);
void SphericalLatlonToCartesian(Real *x, Real *y, Real *z, Real a, Real b, Real c, Real phi, Real theta);
void CartesianToSphericalLatlon(Real *a, Real *b, Real *c, Real x, Real y, Real z, Real phi, Real theta);

// History output total absolute angular momentum and energy
Real TotalAbsoluteAngularMomentum(MeshBlock *pm, int iout);
Real TotalEnergy(MeshBlock *pm, int iout);

// particle functions
bool ParticleTranslate(MeshBlock *pmb, Particle &pt, int cid[3], Real const time, Real const dt);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalAbsoluteAngularMomentum, "AM");
  EnrollUserHistoryOutput(1, TotalEnergy, "EN");
  EnrollUserParticleUpdateFunction(ParticleTranslate);
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // read and save problem parameter
  int  vnum = pin->GetInteger("problem", "vnum");
  Real vrad = pin->GetReal("problem", "vrad");
  Real vlat = pin->GetReal("problem", "vlat")/180.*M_PI;
  Real vgh  = pin->GetReal("problem", "vgh");
  Real gh0  = pin->GetReal("problem", "gh0");
  int  ntracers = pin->GetInteger("problem", "ntracers");

  Real radius = pcoord->x3v(0);

  // setup vortex longitude
  Real *vlon = new Real [vnum];
  for (int n = 0; n < vnum; ++n)
    vlon[n] = 2.*M_PI*n/vnum;

  // setup initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real theta, phi, xx, yy, zz, u1, u2, u3;
        RotateEquatorToPole(&theta, &phi, pcoord->x2v(j), pcoord->x1v(i));

        // add vortices around the pole
        // vortex center position is C = [vlat, 2*pi*n/vnum]
        // vortex height is gh_vortex * exp(-0.5*sqr(X - C)/sqr(vrad))
        phydro->u(IDN,k,j,i) = gh0;
        for (int n = 0; n < vnum; ++n) {
          Real dist = radius * acos(sin(vlat)*sin(theta) + cos(vlat)*cos(theta)*cos(vlon[n]-phi));
          phydro->u(IDN,k,j,i) += vgh * exp(-0.5*_sqr(dist)/_sqr(vrad));
        }

        // another vortex at the pole
        Real dist = radius * (M_PI/2. - theta);
        phydro->u(IDN,k,j,i) += vgh * exp(-0.5*_sqr(dist)/_sqr(vrad));
      }

  // randomly distribute tracers over the domain
  ppg = new ParticleGroup(this, "tracer");
  long int iseed = -1 - Globals::my_rank;

  Real x1min = block_size.x1min*0.99;
  Real x1max = block_size.x1max*0.99;
  Real x2min = block_size.x2min*0.99;
  Real x2max = block_size.x2max*0.99;

  /*for (int n = 0; n < ntracers; ++n) {
    Particle tracer;
    tracer.time = 0.;
    tracer.x3 = 0.;
    tracer.x2 = asin(sin(x2min) + (sin(x2max) - sin(x2min))*ran2(&iseed));
    tracer.x1 = x1min + (x1max - x1min) * ran2(&iseed);
    ppg->q.push_back(tracer);
  }*/

  Real lat, lon;
  for (int n2 = 0; n2 < sqrt(ntracers); ++n2)
    for (int n1 = 0; n1 < sqrt(ntracers); ++n1) {
      Particle tracer;
      tracer.x1 = x1min + (x1max - x1min)*(n1/sqrt(ntracers));
      tracer.x2 = asin(sin(x2min) + (sin(x2max) - sin(x2min))*(n2/sqrt(ntracers)));
      tracer.x3 = 0.;
      RotateEquatorToPole(&lat, &lon, tracer.x2, tracer.x1);
      tracer.time = lat;
      ppg->q.push_back(tracer);
  }

  delete[] vlon;
}

bool ParticleTranslate(MeshBlock *pmb, Particle &pt, int cid[3], Real const time, Real const dt)
{
  Real radius = pmb->pcoord->x3v(0),
       lat = pmb->pcoord->x2v(cid[1]);
  pt.x1 += pt.v1 * dt / (radius * cos(lat));
  pt.x2 += pt.v2 * dt / radius;
  return true;
}

void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0)
{
  *theta1 = asin(cos(theta0)*sin(phi0));
  *phi1   = asin(cos(theta0)*cos(phi0)/cos(*theta1));

  if (theta0 < 0. && (phi0 > -M_PI/2. && phi0 < M_PI/2.))
    *phi1 = M_PI - *phi1;
  if (theta0 < 0. && (phi0 < -M_PI/2. || phi0 > M_PI/2.))
    *phi1 = - M_PI - *phi1;
}

void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0)
{
  *theta1 = asin(cos(theta0)*cos(phi0));
  *phi1   = asin(sin(theta0)/cos(*theta1));

  if (phi0 < 0. && theta0 > 0.)
    *phi1 = M_PI - *phi1;
  if (phi0 < 0. && theta0 < 0.)
    *phi1 = - M_PI - *phi1;
}

void SphericalLatlonToCartesian(
    Real *x, Real *y, Real *z, 
    Real a, Real b, Real c,
    Real phi, Real theta)
{
  *x = -a*sin(phi) - b*sin(theta)*cos(phi) + c*cos(theta)*cos(phi);
  *y = a*cos(phi) - b*sin(theta)*sin(phi) + c*cos(theta)*sin(phi);
  *z = b*cos(theta) + c*sin(theta);
}

void CartesianToSphericalLatlon(
    Real *a, Real *b, Real *c, 
    Real x, Real y, Real z, 
    Real phi, Real theta)
{
  *a = -x*sin(phi) + y*cos(phi);
  *b = -x*sin(theta)*cos(phi) - y*sin(theta)*sin(phi) + z*cos(theta);
  *c = x*cos(theta)*cos(phi) + y*cos(theta)*sin(phi) + z*sin(theta);
}

Real TotalAbsoluteAngularMomentum(MeshBlock *pmb, int iout)
{
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);
  Hydro *phyd = pmb->phydro;
  Real am = 0.;
  Real omega = phyd->psrc->GetOmegaZ();

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real phi = phyd->u(IDN,k,j,i);
        Real r = pmb->pcoord->x3v(k);
        Real theta = pmb->pcoord->x2v(j);
        Real u = phyd->u(IM1,k,j,i);
        am += vol(i) * (omega*phi*r*cos(theta) + u) * r*cos(theta);
      }
    }

  return am;
}

Real TotalEnergy(MeshBlock *pmb, int iout)
{
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);
  Hydro *phyd = pmb->phydro;
  Real en = 0.;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real phi = phyd->u(IDN,k,j,i);
        Real ke1 = 0.5 * _sqr(phyd->u(IM1,k,j,i)) / phi;
        Real ke2 = 0.5 * _sqr(phyd->u(IM2,k,j,i)) / phi;
        Real pe = 0.5 * _sqr(phi);
        en += vol(i) * (ke1 + ke2 + pe);
      }
    }

  return en;
}
