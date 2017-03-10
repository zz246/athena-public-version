//! \file  roe_shallow_water.cpp
//  \brief Roe's linearized Riemann solver for shallow water model

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>

// Athena++ headers
#include "../hydro.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  int ivy = ivx == 1 ? 2 : 1;

  Real ul[3], ur[3], wave[3][3], speed[3];

  Real ubar, vbar, cbar, delh, deluh, delvh, a1, a2, a3;

  for (int i = il; i <= iu; ++i) {
    ul[0] = wl(IDN, i);
    ul[1] = wl(IDN, i) * wl(ivx, i);
    ul[2] = wl(IDN, i) * wl(ivy, i);

    ur[0] = wr(IDN, i);
    ur[1] = wr(IDN, i) * wr(ivx, i);
    ur[2] = wr(IDN, i) * wr(ivy, i);

    ubar = (ul[1] / sqrt(ul[0]) + ur[1] / sqrt(ur[0])) / (sqrt(ul[0]) + sqrt(ur[0]));
    vbar = (ul[2] / sqrt(ul[0]) + ur[2] / sqrt(ur[0])) / (sqrt(ul[0]) + sqrt(ur[0]));
    cbar = sqrt(0.5 * (ul[0] + ur[0]));

    delh = ur[0] - ul[0];
    deluh = ur[1] - ul[1];
    delvh = ur[2] - ul[2];

    a1 = 0.5 * (-deluh + (ubar + cbar) * delh) / cbar;
    a2 = delvh - vbar * delh;
    a3 = 0.5 * ( deluh - (ubar - cbar) * delh) / cbar;

    wave[0][0] = a1; wave[0][1] = a1 * (ubar - cbar); wave[0][2] = a1 * vbar;
    wave[1][0] = 0.; wave[1][1] = 0.; wave[1][2] = a2;
    wave[2][0] = a3; wave[2][1] = a3 * (ubar + cbar); wave[2][2] = a3 * vbar;

    speed[0] = ubar - cbar;
    speed[1] = ubar;
    speed[2] = ubar + cbar;

    flx(IDN, i) = ul[1];
    flx(ivx, i) = ul[1] * ul[1] / ul[0] + 0.5 * ul[0] * ul[0];
    flx(ivy, i) = ul[1] * ul[2] / ul[0];
    flx(IVZ, i) = 0.;

    bool entropy_fix = false;

    if (!entropy_fix) {
      for (int r = 0; r < 3; ++r)
        if (speed[r] < 0.) {
          flx(IDN, i) += speed[r] * wave[r][0];
          flx(ivx, i) += speed[r] * wave[r][1];
          flx(ivy, i) += speed[r] * wave[r][2];
        }
    } else {
      Real s0, s1, sfract, tmp[3];

      // check 1-wave
      s0 = ul[1]/ul[0] - sqrt(ul[0]);
      for (int r = 0; i < 3; ++r)
        tmp[r] = ul[r] + wave[0][r];
      s1 = tmp[1]/tmp[0] - sqrt(tmp[0]);
      if (s0 < 0. && s1 > 0.)
        sfract = s0 * (s1 - speed[0])/(s1 - s0);
      else if (speed[0] < 0.)
        sfract = speed[0];
      else 
        sfract = 0.;

      flx(IDN, i) += sfract * wave[0][0];
      flx(ivx, i) += sfract * wave[0][1];
      flx(ivy, i) += sfract * wave[0][2];

      // check 2-wave
      if (speed[1] < 0.)
        sfract = speed[1];
      else
        sfract = 0.;

      flx(IDN, i) += sfract * wave[1][0];
      flx(ivx, i) += sfract * wave[1][1];
      flx(ivy, i) += sfract * wave[1][2];

      // check 3-wave
      s1 = ur[1]/ur[0] + sqrt(ur[0]);
      for (int r = 0; r < 3; ++r)
        tmp[r] = ur[r] - wave[2][r];
      s0 = tmp[1]/tmp[0] + sqrt(tmp[0]);

      if (s0 < 0. && s1 > 0.)
        sfract = s0 * (s1 - speed[2])/(s1 - s0);
      else if (speed[2] < 0.)
        sfract = speed[2];
      else
        sfract = 0.;

      flx(IDN, i) += sfract * wave[2][0];
      flx(ivx, i) += sfract * wave[2][1];
      flx(ivy, i) += sfract * wave[2][2];
    }
  }
}
