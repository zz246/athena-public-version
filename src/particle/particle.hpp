#ifndef PARTICLE_HPP
#define PARTICLE_HPP

// C++ headers
#include <vector>
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
class ParticleTableOutput;

struct Particle {
  Real time, x1, x2, x3;
  Real v1, v2, v3;

  #if NREAL_PARTICLE_DATA > 0
    Real rdata[NREAL_PARTICLE_DATA];
  #endif

  #if NINT_PARTICLE_DATA > 0
    int  idata[NINT_PARTICLE_DATA];
  #endif
};

#ifdef MPI_PARALLEL
extern MPI_Datatype MPI_PARTICLE;
#endif

class ParticleGroup {
  //friend class ParticleTableOutput;
public:
  ParticleGroup(MeshBlock *pmb, std::string _name);
  ~ParticleGroup();
  
  // data
  MeshBlock* pmy_block;
  std::string name;
  ParticleGroup *prev, *next;
  std::vector<Particle> q;
  std::vector<int> bufid;
  static int ntotal;

  // functions
  ParticleGroup* AddParticleGroup(MeshBlock *pmb, std::string name);
  std::vector<Particle>& GetParticle(std::string name);
  std::vector<Particle> const& GetParticle(std::string name) const;
  void PropertyUpdate(Real time, Real dt);

protected:
  ParticleUpdateFunc_t particle_fn_;
  Real *coordinates_;
  int lengths_[3];
};

#endif
