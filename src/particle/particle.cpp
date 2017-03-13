// Athena++ headers
#include "particle.hpp"

// constructor, initializes data structure and parameters
Particle::Particle(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  prev = NULL;
  next = NULL;
}

// destructor
Particle::~Particle()
{}

// functions
void InterpolateToCellCenter(AthenaArray<Real> &u, Coordinates *pcoord)
{}

void InterpolateToCellFace(AthenaArray<Real> &u, Coordinates *pcoord)
{}
