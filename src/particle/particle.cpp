// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "particle.hpp"

// constructor, initializes data structure and parameters
ParticleGroup::ParticleGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name)
{
  prev = NULL;
  next = NULL;
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
  ParticleGroup *p = this;
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

  return p->data;
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

  return p->data;
}

void ParticleGroup::InterpolateToCellCenter(AthenaArray<Real> &u, Coordinates *pcoord) const
{}

void ParticleGroup::InterpolateToCellFace(AthenaArray<Real> &u, Coordinates *pcoord) const
{}
