// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "particle.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../athena_math.hpp" // _interpn

// constructor, initializes data structure and parameters
ParticleGroup::ParticleGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name)
{
  prev = NULL;
  next = NULL;
  particle_fn_ = pmb->pmy_mesh->particle_fn_;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  coordinates_.resize(ncells3 + ncells2 + ncells1);

  for (int k = 0; k < ncells3; ++k)
    coordinates_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < ncells2; ++j)
    coordinates_[ncells3 + j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < ncells1; ++i)
    coordinates_[ncells3 + ncells2 + i] = pmb->pcoord->x1v(i);

  lengths_[0] = ncells3;
  lengths_[1] = ncells2;
  lengths_[2] = ncells1;
}

// destructor
ParticleGroup::~ParticleGroup()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

// functions
ParticleGroup* ParticleGroup::AddParticleGroup(MeshBlock *pmb, std::string name)
{
  std::stringstream msg;
  ParticleGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddParticleGroup: ParticleGroup is empty, use new ParticleGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ParticleGroup(pmb, name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

std::vector<Particle>& ParticleGroup::GetParticle(std::string name)
{
  std::stringstream msg;
  ParticleGroup *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}

std::vector<Particle> const& ParticleGroup::GetParticle(std::string name) const
{
  std::stringstream msg;
  ParticleGroup const *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}

void ParticleGroup::PropertyUpdate(Real time, Real dt)
{
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];

  Hydro *phydro = pmy_block->phydro;

  v1.InitWithShallowSlice(phydro->w,4,IM1,1);
  v2.InitWithShallowSlice(phydro->w,4,IM2,1);
  v3.InitWithShallowSlice(phydro->w,4,IM3,1);

  Real x1min = pmy_block->block_size.x1min;
  Real x1max = pmy_block->block_size.x1max;
  Real x2min = pmy_block->block_size.x2min;
  Real x2max = pmy_block->block_size.x2max;
  Real x3min = pmy_block->block_size.x3min;
  Real x3max = pmy_block->block_size.x3max;

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

    particle_fn_(q[i], time, dt);

    int ox1 = q[i].x1 < x1min ? -1 : (q[i].x1 > x1max ? 1 : 0);
    int ox2 = q[i].x2 < x2min ? -1 : (q[i].x2 > x2max ? 1 : 0);
    int ox3 = q[i].x3 < x3min ? -1 : (q[i].x3 > x3max ? 1 : 0);
  }
}
