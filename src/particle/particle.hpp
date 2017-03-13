#ifndef PARTICLE_HPP
#define PARTICLE_HPP

// C++ headers
#include <vector>

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParameterInput;

struct OneParticle {
  Real x1, x2, x3, time;
  Real v1, v2, v3;
  Real rdata[NREAL_PARTICLE_DATA];
  int  idata[NINT_PARTICLE_DATA];
};

class Particle {
public:
  Particle(MeshBlock *pmb, ParameterInput *pin);
  ~Particle();

  // data
  MeshBlock* pmy_block;
  std::vector<OneParticle> q;

  // functions
  void InterpolateToCellCenter(AthenaArray<Real> &u, Coordinates *pcoord);
  void InterpolateToCellFace(AthenaArray<Real> &u, Coordinates *pcoord);

  Particle *prev, *next;
};

#endif
