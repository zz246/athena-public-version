#ifndef AIRPARCEL_HPP_
#define AIRPARCEL_HPP_

// C++ headers
#include <iosfwd>

// Athena++ headers
#include "../particle/particle.hpp"

/** @file
 * @brief This file contains declarations of AirParcel and AirParcelD.
 * 
 * **Author**   : Cheng Li, California Institute of Technology <br>
 * **Contact**  : cli@gps.caltech.edu <br>
 * **Revision history** :
 * - Sep 3 2015, first version.
 * - Oct 27 2015, change set_n_dry_moist to static
 * - Feb 11 2016, add three functions that operates on AirParcel:
 *  + normalize_to_mixr,
 *  + set_air_from,
 *  + set_air_to
 * - April 26 2016, Version 3
 */

/** @brief AirParcel is a light-weighted class that describes an air parcel.
 *
 *  It collects the basic thermodynamic properties of an air parcel, i.e.
 *  temperature, pressure and mixing ratios of each species (ordered by dry, moist 
 *  ,and condensates). It uses an cache-friendly memory layout, in which
 *  mixing ratios are stored in fixed sized array (std::array<double, MAXMOL>). 
 *  MAXMOL is a macro defined in this file that specifies the maximum 
 *  number of species to hold.
 *
 *  The reason to abstract AirParcel out of thermodynamics
 *  calculations is that, in many cases, only several basic properties of the air
 *  parcel are needed. For example, dynamics does not care about the entropy and
 *  enthalpy of an air parcel but only needs the mixing ratios. AirParcel is the 
 *  base class that contains these information.
 */

class AirParcel : protected Particle {
  friend std::ostream& operator<<(std::ostream &os, AirParcel const& air);
public:
  AirParcel();

  // data
  static int ndry, ngas, ntotal;

  // functions
  void Normalize();
  void SetFrom(double const* M, double T, double P);
  void SetTo(double *M, double *T, double *P);
  void SetTP(double temp, double pres) {
      rdata[0] = temp;
      rdata[1] = pres;
  }

  double GetTemp() const { return rdata[0]; }
  double GetPres() const { return rdata[1]; }
  double GetBar() const { return 1.E-5 * rdata[1]; }

  double& Mol(int i) { return rdata[2 + i]; }
  double const& Mol(int i) const { return rdata[2 + i]; }

  AirParcel& operator+=(AirParcel const& air);
  AirParcel& operator-=(AirParcel const& air);
  AirParcel& operator*=(AirParcel const& air);
  AirParcel& operator/=(AirParcel const& air);
  AirParcel& operator*=(double r);
  AirParcel& operator/=(double r);
  AirParcel operator+(AirParcel const& air) const;
  AirParcel operator-(AirParcel const& air) const;
  AirParcel operator*(AirParcel const& air) const;
  AirParcel operator/(AirParcel const& air) const;
  AirParcel operator*(double r) const;
  AirParcel operator/(double r) const;
    
  double Mall(int i) const;
};

AirParcel operator*(double r, AirParcel const& air);

#endif
