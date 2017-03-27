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
#include "../utils/utils.hpp"
#include "../particle/particle.hpp"

// model parameters shared by all subroutine
namespace Prob {
  static Real r_storm;
  static Real tau_interval;
  static Real tau_storm;
  static Real smax;
  static Real tau_mass;
  static Real tau_ape;
  static Real gheq;
  static int  tracers_per_storm;
  static Real tracer_lifetime;

  static Real theta1;
  static Real theta2;
  static Real umax;

  static Real radius;
  static Real omega;
}

// last time of mass pulse in this MeshBlock
static Real t_last_storm = -1.E-8;

// seed of random number
static long int iseed = -1;

// external forcing
void MassPulseNewtonianCooling(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void LagrangianTracerAdvection(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// support functions
Real GetZonalWind(Real theta);
Real SolveFarFieldGeopotential(Real gh0);
Real GetGeopotential(Real gh0, Real theta);
void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0);
void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0);
void SphericalLatlonToCartesian(Real *x, Real *y, Real *z, Real a, Real b, Real c, Real phi, Real theta);
void CartesianToSphericalLatlon(Real *a, Real *b, Real *c, Real x, Real y, Real z, Real phi, Real theta);
void SetTracerLatlon(Real *lat, Real *lon, Real theta, Real phi, Real alpha, Real delta);

// History output total absolute angular momentum and energy
Real TotalAbsoluteAngularMomentum(MeshBlock *pm, int iout);
Real TotalEnergy(MeshBlock *pm, int iout);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalAbsoluteAngularMomentum, "AM");
  EnrollUserHistoryOutput(1, TotalEnergy, "EN");

  //EnrollUserExplicitSourceFunction(MassPulseNewtonianCooling);
  EnrollUserExplicitSourceFunction(LagrangianTracerAdvection);
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // read and save problem parameter
  Prob::r_storm       = pin->GetReal("problem", "r_storm");
  Prob::tau_interval  = pin->GetReal("problem", "tau_interval");
  Prob::tau_storm     = pin->GetReal("problem", "tau_storm");
  Prob::smax          = pin->GetReal("problem", "smax");
  Prob::tau_mass      = pin->GetReal("problem", "tau_mass");
  Prob::tau_ape       = pin->GetReal("problem", "tau_ape");
  Prob::gheq          = pin->GetReal("problem", "gheq");
  Prob::tracers_per_storm = pin->GetInteger("problem", "tracers_per_storm");
  Prob::tracer_lifetime = pin->GetReal("problem", "tracer_lifetime");

  Prob::theta1        = pin->GetReal("problem", "theta1");
  Prob::theta2        = pin->GetReal("problem", "theta2");
  Prob::umax          = pin->GetReal("problem", "umax");

  Prob::radius        = pcoord->x3v(0);
  Prob::omega         = phydro->psrc->GetOmegaX();

  // solve for far field geopotential
  Real gh0;
  int err = _root(Prob::gheq/2., Prob::gheq*2., 1., &gh0, SolveFarFieldGeopotential);

  if (err != 0)
    throw std::runtime_error("FATAL ERROR: SolveFarFieldGeopotential does not converge");

  // setup initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real theta, phi, xx, yy, zz, u1, u2, u3;
        RotateEquatorToPole(&theta, &phi, pcoord->x2v(j), pcoord->x1v(i));

        // setup mean flow
        phydro->u(IDN,k,j,i) = GetGeopotential(gh0, theta);

        SphericalLatlonToCartesian(&xx,&yy,&zz,GetZonalWind(theta),0.,0.,phi,theta);
        CartesianToSphericalLatlon(&u1,&u2,&u3,zz,xx,yy,pcoord->x1v(i),pcoord->x2v(j));
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i) * u1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i) * u2;

        // save initial geopotential field
        //gh_init(k,j,i) = 

        // add pertubation
        phydro->u(IDN,k,j,i) += 2400.*cos(theta)*exp(-_sqr(3.*phi))
          * exp(-_sqr(15.*((Prob::theta1+Prob::theta2)/2.-theta)));
      }

  // storm and tracer particles
  //ppg = new ParticleGroup(this, "storm");
  //ppg->AddParticleGroup(this, "tracer");
  
  // insert randomly distributed tracers over the domain
  ppg = new ParticleGroup(this, "tracer");

  Real lonmin = pmy_mesh->mesh_size.x1min;
  Real lonmax = pmy_mesh->mesh_size.x1max;
  Real latmin = pmy_mesh->mesh_size.x2min;
  Real latmax = pmy_mesh->mesh_size.x2max;

  std::vector<Particle>& tracers = ppg->GetParticle("tracer");
  for (int n = 0; n < Prob::tracers_per_storm; ++n) {
    Particle tracer;
    tracer.time = 0.;
    tracer.x2 = asin(sin(latmin) + (sin(latmax) - sin(latmin))*ran2(&iseed));
    tracer.x1 = lonmin + (lonmax - lonmin)*ran2(&iseed);
    tracers.push_back(tracer);
  }
}

void LagrangianTracerAdvection(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real lonmin = pmb->pmy_mesh->mesh_size.x1min;
  Real lonmax = pmb->pmy_mesh->mesh_size.x1max;
  Real latmin = pmb->pmy_mesh->mesh_size.x2min;
  Real latmax = pmb->pmy_mesh->mesh_size.x2max;

  // Lagrangian tracer advection and elimination
  std::vector<Particle>& tracers = pmb->ppg->GetParticle("tracer");

  if (step == 2) {
    pmb->phydro->InterpolateVelocity(tracers);
    for (size_t n = 0; n < tracers.size(); ++n) {
      Real age = time - tracers[n].time;
      if (tracers[n].x2 < latmin || tracers[n].x2 > latmax ||
          tracers[n].x1 < lonmin || tracers[n].x1 > lonmax ||
          ran2(&iseed) < exp(-Prob::tracer_lifetime / age)) {
        tracers[n].time = -1.;
      } else {
        tracers[n].x1 += tracers[n].v1 * dt / (Prob::radius * cos(tracers[n].x2));
        tracers[n].x2 += tracers[n].v2 * dt / Prob::radius;
      }
    }

    size_t n = 0;
    while (n < tracers.size()) {
      if (tracers[n].time < 0.)
        tracers.erase(tracers.begin() + n);
      else
        n++;
    }
  }
}

void MassPulseNewtonianCooling(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  std::vector<Particle>& storms = pmb->ppg->GetParticle("storm");
  std::vector<Particle>& tracers = pmb->ppg->GetParticle("tracer");

  size_t current_storms = storms.size();
  Real lonmin = pmb->pmy_mesh->mesh_size.x1min;
  Real lonmax = pmb->pmy_mesh->mesh_size.x1max;
  Real latmin = pmb->pmy_mesh->mesh_size.x2min;
  Real latmax = pmb->pmy_mesh->mesh_size.x2max;

  // insert new storm with tracer
  if ((step == 1) && (ran2(&iseed) < exp(-Prob::tau_interval/(time - t_last_storm)))) {
    Particle storm;
    storm.time = time;
    storm.x2 = asin(sin(latmin) + (sin(latmax) - sin(latmin))*ran2(&iseed));
    storm.x1 = lonmin + (lonmax - lonmin)*ran2(&iseed);
    storms.push_back(storm);
    t_last_storm = time;

    for (int n = 0; n < Prob::tracers_per_storm; ++n) {
      Particle tracer;
      Real alpha = asin(2.2*Prob::r_storm/Prob::radius)*ran2(&iseed);
      Real delta = 2.*M_PI*ran2(&iseed);

      tracer.time = time;
      SetTracerLatlon(&tracer.x2, &tracer.x1, storm.x2, storm.x1, alpha, delta);
      tracers.push_back(tracer);
    }
  }

  /* do not use newtonian cooling now
  // calculate average height for this instant
  Real gh_avg[2];
  AthenaArray<Real> vol;
  int ncells1 = pmb->block_size.nx1 + 2*NGHOST;
  vol.NewAthenaArray(ncells1);

  gh_avg[0] = 0.;
  gh_avg[1] = 0.;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        gh_avg[0] += prim(IDN,k,j,i) * vol(i);
        gh_avg[1] += vol(i);
      }
    }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, gh_avg, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  gh_avg[0] /= gh_avg[1];
  */

  // add storm source
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lon1 = pmb->pcoord->x1v(i);
        Real lat1 = pmb->pcoord->x2v(j);
        for (size_t n = 0; n < current_storms; ++n) {
          Real lon2 = storms[n].x1;
          Real lat2 = storms[n].x2;
          Real dr = Prob::radius * acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
          if (dr > 2.2 * Prob::r_storm) continue;

          // mass pulse
          Real r2r2 = _sqr(dr)/_sqr(Prob::r_storm);
          Real t2t2 = _sqr(time - (storms[n].time + Prob::tau_storm/2.))/_sqr(Prob::tau_storm);
          cons(IDN,k,j,i) += dt * Prob::smax * exp(- r2r2 - t2t2);
        }
        // Newtonian cooling
        //cons(IDN,k,j,i) += - dt * (gh_avg[0] - gheq)/tau_mass;
      }

  // remove deceased storms
  if ((step == 2) && (storms.size() > 0)) {
    size_t n = 0;
    for (; n < storms.size(); ++n)
      if (time - storms[n].time < Prob::tau_storm)
        break;
    storms.erase(storms.begin(), storms.begin() + n);
  }

  // Lagrangian tracer advection and elimination
  if (step == 2) {
    pmb->phydro->InterpolateVelocity(tracers);
    for (size_t n = 0; n < tracers.size(); ++n) {
      Real age = time - tracers[n].time;
      if (tracers[n].x2 < latmin || tracers[n].x2 > latmax ||
          tracers[n].x1 < lonmin || tracers[n].x1 > lonmax ||
          ran2(&iseed) < exp(-Prob::tracer_lifetime / age)) {
        tracers[n].time = -1.;
      } else {
        tracers[n].x1 += tracers[n].v1 * dt / (Prob::radius * cos(tracers[n].x2));
        tracers[n].x2 += tracers[n].v2 * dt / Prob::radius;
      }
    }

    size_t n = 0;
    while (n < tracers.size()) {
      if (tracers[n].time < 0.)
        tracers.erase(tracers.begin() + n);
      else
        n++;
    }
  }
}

Real GetZonalWind(Real theta)
{
  if (theta <= Prob::theta1 || theta >= Prob::theta2)
    return 0.;

  Real en = exp(-4./_sqr(Prob::theta2 - Prob::theta1));

  return Prob::umax/en * exp(1./((theta - Prob::theta1) * (theta - Prob::theta2)));
}

Real GetGeopotential(Real gh0, Real theta)
{
  Real ds = 0.001;
  for (Real s = - M_PI/2. + ds; s < theta; s += ds) {
    Real f1r = 2.*Prob::omega*sin(s-ds)*Prob::radius,
         f2r = 2.*Prob::omega*sin(s)*Prob::radius;

    Real u1 = GetZonalWind(s-ds),
         u2 = GetZonalWind(s);

    gh0 -= 0.5 * ((f1r*u1 + u1*u1*tan(s-ds)) + (f2r*u2 + u2*u2*tan(s))) * ds;
  }

  return gh0;
}

Real SolveFarFieldGeopotential(Real gh)
{
  Real total_gh = 0.;
  Real ds = 0.01, theta = - M_PI/2. + ds;

  for (; theta < M_PI/2.; theta += ds) {
    Real f1r = 2.*Prob::omega*sin(theta-ds)*Prob::radius,
         f2r = 2.*Prob::omega*sin(theta)*Prob::radius;

    Real u1 = GetZonalWind(theta-ds),
         u2 = GetZonalWind(theta);

    gh -= 0.5 * ((f1r*u1 + u1*u1*tan(theta-ds)) + (f2r*u2 + u2*u2*tan(theta))) * ds;
    total_gh += gh * cos(theta) * ds;
  }

  return total_gh/2. - Prob::gheq;
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

void SetTracerLatlon(Real *lat, Real *lon, Real theta, Real phi, Real alpha, Real delta)
{
  *lat = asin(sin(theta)*cos(alpha) - cos(theta)*sin(alpha)*cos(delta));
  *lon = phi + acos((sin(theta)*sin(alpha)*cos(delta) + cos(theta)*cos(alpha))/cos(*lat));
  if (*lon > M_PI) *lon -= 2.*M_PI;
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
