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

// model parameters shared by all subroutine
namespace Prob {
  static Real r_storm;
  static Real tau_interval;
  static Real tau_storm;
  static Real smax;
  static Real tau_mass;
  static Real tau_ape;
  static Real gheq;

  static Real theta1;
  static Real theta2;
  static Real umax;

  static Real radius;
  static Real omega;
}


// mass pulse
struct MassPulse {
  Real tpeak;
  Real lat;
  Real lon;
  int  cyclic;
};

// all mass pulses
static std::vector<MassPulse> storms;

// last time of mass pulse in this MeshBlock
static Real t_last_storm = -1.E-8;

// seed of random number
static long int iseed = -1;

// external forcing
void MassPulseRadiationSink(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// support functions
Real GetZonalWind(Real theta);
Real SolveFarFieldGeopotential(Real gh0);
Real GetGeopotential(Real gh0, Real theta);
void RotatePoleToEquator(Real *theta1, Real *phi1, Real theta0, Real phi0);
void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0);
void SphericalLatlonToCartesian(Real *x, Real *y, Real *z, Real a, Real b, Real c, Real phi, Real theta);
void CartesianToSphericalLatlon(Real *a, Real *b, Real *c, Real x, Real y, Real z, Real phi, Real theta);

// History output total absolute angular momentum and energy
Real TotalAbsoluteAngularMomentum(MeshBlock *pm, int iout);
Real TotalEnergy(MeshBlock *pm, int iout);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, TotalAbsoluteAngularMomentum, "AM");
  EnrollUserHistoryOutput(1, TotalEnergy, "EN");

  //EnrollUserExplicitSourceFunction(MassPulseRadiationSink);
}

//  \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Prob::r_storm       = pin->GetReal("problem", "r_storm");
  Prob::tau_interval  = pin->GetReal("problem", "tau_interval");
  Prob::tau_storm     = pin->GetReal("problem", "tau_storm");
  Prob::smax          = pin->GetReal("problem", "smax");
  Prob::tau_mass      = pin->GetReal("problem", "tau_mass");
  Prob::tau_ape       = pin->GetReal("problem", "tau_ape");
  Prob::gheq          = pin->GetReal("problem", "gheq");

  Prob::theta1        = pin->GetReal("problem", "theta1");
  Prob::theta2        = pin->GetReal("problem", "theta2");
  Prob::umax          = pin->GetReal("problem", "umax");

  Prob::radius        = pcoord->x3v(0);
  Prob::omega         = phydro->psrc->GetOmegaX();

  // solve for far field geopotential
  Real gh0;
  int err = _root(Prob::gheq/2., Prob::gheq*2., 1., &gh0, SolveFarFieldGeopotential);
  std::cout << gh0 << std::endl;

  if (err != 0)
    throw std::runtime_error("FATAL ERROR: SolveFarFieldGeopotential does not converge");

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

        // add pertubation
        phydro->u(IDN,k,j,i) += 2400.*cos(theta)*exp(-_sqr(3.*phi))
          * exp(-_sqr(15.*((Prob::theta1+Prob::theta2)/2.-theta)));
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
}

void RotateEquatorToPole(Real *theta1, Real *phi1, Real theta0, Real phi0)
{
  *theta1 = asin(cos(theta0)*cos(phi0));
  *phi1   = asin(sin(theta0)/cos(*theta1));

  if (phi0 < 0. && theta0 > 0.)
    *phi1 = M_PI - *phi1;
  if (phi0 < 0. && theta0 < 0.)
    *phi1 = - *phi1 - M_PI;
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

/*void MassPulseRadiationSink(MeshBlock *pmb, const Real time, const Real dt, const int step,
      const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  int current_storms = storms.size();
  Real lonmin = pmb->pmy_mesh->mesh_size.x1min;
  Real lonmax = pmb->pmy_mesh->mesh_size.x1max;
  Real latmin = pmb->pmy_mesh->mesh_size.x2min;
  Real latmax = pmb->pmy_mesh->mesh_size.x2max;

  // insert new storm
  if ((step == 1) && (ran2(&iseed) < exp(-tau_interval/(time - t_last_storm)))) {
    MassPulse storm;
    storm.tpeak = time + tau_storm/2.;
    storm.lat = asin(sin(latmin) + (sin(latmax) - sin(latmin))*ran2(&iseed));
    storm.lon = lonmin + (lonmax - lonmin)*ran2(&iseed);
    storm.cyclic = ran2(&iseed) > 0.5 ? 1 : -1;
    storms.push_back(storm);
    t_last_storm = time;
  }

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

  // add storm source
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real lon1 = pmb->pcoord->x1v(i);
        Real lat1 = pmb->pcoord->x2v(j);
        for (size_t n = 0; n < current_storms; ++n) {
          Real lon2 = storms[n].lon;
          Real lat2 = storms[n].lat;
          Real dr = Prob::radius*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
          int  sign = storms[n].cyclic;
          if (dr > 2.2 * r_storm) continue;

          // mass pulse
          cons(IDN,k,j,i) += dt * smax * sign * exp(-_sqr(dr)/_sqr(r_storm) - _sqr(time - storms[n].tpeak)/_sqr(tau_storm));
        }
        // radiation
        //cons(IDN,k,j,i) += - dt * ((gh_avg[0] - gheq)/tau_mass + (prim(IDN,k,j,i) - gh_avg[0])/tau_ape);
        cons(IDN,k,j,i) += - dt * (gh_avg[0] - gheq)/tau_mass;
      }

  // remove deceased storms
  if ((step == 2) && (storms.size() > 0)) {
    size_t n = 0;
    for (; n < storms.size(); ++n)
      if (time - storms[n].tpeak < 2.2*tau_storm)
        break;
    storms.erase(storms.begin(), storms.begin() + n);
  }
}*/
