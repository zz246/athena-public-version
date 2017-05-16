// C++ headers
#include <iostream>

// Athena++ headers
#include "air_parcel.hpp"

int AirParcel::ndry = 0;
int AirParcel::ngas = 0;
int AirParcel::ntotal = 0;

std::ostream& operator<<(std::ostream &os, AirParcel const& air)
{
  os << "T = " << air.rdata[0] << ", P = " << air.rdata[1] << std::endl;
  for (int i = 0; i < air.ntotal; ++i)
    os << air.Mol(i) << " ";
  os << std::endl;
  return os;
}

AirParcel::AirParcel()
{
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    rdata[i] = 0.;
}

void AirParcel::Normalize() {
  double mols = 0;
  for (int i = 0; i < ntotal; ++i) mols += Mol(i);
  for (int i = 0; i < ntotal; ++i) Mol(i) /= mols;
}

void AirParcel::SetFrom(double const* M, double T, double P) {
  for (int i = 0; i < ntotal; ++i)
    Mol(i) = M[i];
  SetTP(T, P);
}

void AirParcel::SetTo(double *M, double *T, double *P) {
  for (int i = 0; i < ntotal; ++i)
    M[i] = Mol(i);
  *T = GetTemp();
  *P = GetPres();
}

AirParcel& AirParcel::operator+=(AirParcel const& air)
{
  rdata[0] += air.rdata[0];
  rdata[1] += air.rdata[1];
  for (int i = 0; i < ntotal; ++i)
    Mol(i) += air.Mol(i);
  return *this;
}

AirParcel& AirParcel::operator-=(AirParcel const& air)
{
  rdata[0] -= air.rdata[0];
  rdata[1] -= air.rdata[1];
  for (int i = 0; i < ntotal; ++i)
    Mol(i) -= air.Mol(i);
  return *this;
}

AirParcel& AirParcel::operator*=(AirParcel const& air)
{
  rdata[0] *= air.rdata[0];
  rdata[1] *= air.rdata[1];
  for (int i = 0; i < ntotal; ++i)
    Mol(i) *= air.Mol(i);
  return *this;
}

AirParcel& AirParcel::operator/=(AirParcel const& air)
{
  rdata[0] /= air.rdata[0];
  rdata[1] /= air.rdata[1];
  for (int i = 0; i < ntotal; ++i)
    Mol(i) /= air.Mol(i);
  return *this;
}

AirParcel& AirParcel::operator*=(double r)
{
  rdata[0] *= r;
  rdata[1] *= r;
  for (int i = 0; i < ntotal; ++i)
    Mol(i) *= r;
  return *this;
}

AirParcel& AirParcel::operator/=(double r)
{
  rdata[0] /= r;
  rdata[0] /= r;
  for (int i = 0; i < ntotal; ++i)
    Mol(i) /= r;
  return *this;
}

AirParcel AirParcel::operator+(AirParcel const& air) const
{
  AirParcel other(*this);
  other += air;
 return other;
}

AirParcel AirParcel::operator-(AirParcel const& air) const
{
  AirParcel other(*this);
  other -= air;
  return other;
}

AirParcel AirParcel::operator*(AirParcel const& air) const
{
  AirParcel other(*this);
  other *= air;
  return other;
}

AirParcel AirParcel::operator/(AirParcel const& air) const
{
  AirParcel other(*this);
  other /= air;
  return other;
}

AirParcel AirParcel::operator*(double r) const
{
  AirParcel other(*this);
  other *= r;
  return other;
}

AirParcel AirParcel::operator/(double r) const
{
  AirParcel other(*this);
  other /= r;
  return other;
}

double __attribute__((weak)) AirParcel::Mall(int i) const
{
  return Mol(i);
}

AirParcel operator*(double r, AirParcel const& air)
{
  AirParcel other(air);
  return other *= r;
}
