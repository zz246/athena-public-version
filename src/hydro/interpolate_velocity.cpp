//! \file interpolate_velocity.cpp
//  \brief interpolate mesh velocity to particle velocity

// Athena++
#include "hydro.hpp"
#include "../particle/particle.hpp"
#include "../athena_math.hpp"

void Hydro::InterpolateVelocity(std::vector<OneParticle> &q)
{
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];

  v1.InitWithShallowSlice(w1,4,IM1,1);
  v2.InitWithShallowSlice(w1,4,IM2,1);
  v3.InitWithShallowSlice(w1,4,IM3,1);

  for (size_t i = 0; i < q.size(); ++i) {
    loc[0] = q[i].x3;
    loc[1] = q[i].x2;
    loc[2] = q[i].x1;
    _interpn(&q[i].v1, loc, v1.data(), coordinates_.data(), lengths_, 3);

    if (lengths_[1] > 1)
      _interpn(&q[i].v2, loc, v2.data(), coordinates_.data(), lengths_, 3);
    else
      q[i].v2 = 0.;

    if (lengths_[0] > 1)
      _interpn(&q[i].v3, loc, v3.data(), coordinates_.data(), lengths_, 3);
    else
      q[i].v3 = 0.;
  }
}
