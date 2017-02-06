#ifndef ATHENA_MATH_HPP
#define ATHENA_MATH_HPP
#include <cmath>      // sqrt, fabs
#include <iostream>   // std::endl
#include <sstream>    // std::stringstream
#include <stdexcept>  // std::runtime_error

#include "athena.hpp"  // Real

// small & simple functions
inline Real _sqr(Real x) { return x * x; }
inline Real _max(Real a, Real b) { return a > b ? a : b; }
inline Real _min(Real a, Real b) { return a < b ? a : b; }
inline int _sign(Real x) { return x < 0. ? -1 : 1; }

#undef  MAX_IT
#define MAX_IT 100

#undef  UNLIKELY_VAL
#define UNLIKELY_VAL -1.11111e+30

template<typename T>
int _root(
        Real  x1,
        Real  x2,
        Real  xacc,
        Real *x_root,
        T& func
        )
{
  int
    iter,
    compare;
  Real
    fh,fl,fm,fnew,
    s,xh,xl,xm,xnew;
  static char
    name[]="_root";

  std::stringstream msg;

  fl = func(x1);
  fh = func(x2);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl      = x1;
    xh      = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5*(xl+xh);
      fm = func(xm);
      s  = sqrt(fm*fm-fl*fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm+(xm-xl)*((fl > fh ? 1. : -1.)*fm/s);

      if (fabs(xnew-*x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew    = func(*x_root);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      }
      else {
        msg << "### FATAL ERROR in _root function: should never get here" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if (fabs(xh-xl) <= xacc) {
        return 0;
      }
    }
    msg << "### FATAL ERROR in _root function: exceeded MAX_IT = ";
    msg << MAX_IT;
    msg << "current root calc = ";
    msg << *x_root;
    msg << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  else {
    if (fl == 0.) {
      *x_root = x1;
      return 0;
    }
    if (fh == 0.) {
      *x_root = x2;
      return 0;
    }

    compare = fabs(fl) >= fabs(fh) ? 1 : -1;
    if (compare < 0) {
      return -1;
    }
    else {
      return 1;
    }
  }

  /* Should never get here. */
  msg << "### FATAL ERROR in _root function: should never get here" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

#undef  MAX_IT
#undef  UNLIKELY_VAL

#endif
