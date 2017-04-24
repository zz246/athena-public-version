// Athena++ headers
#include "particle.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../athena_math.hpp" // _interpn

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

  size_t i = 0, j = q.size();
  int ox1, ox2, ox3, fi1, fi2;
  Particle tmp;

  while (i < j) {
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

    int alive = particle_fn_(q[i], time, dt);

    ox1 = q[i].x1 < x1min ? -1 : (q[i].x1 > x1max ? 1 : 0);
    ox2 = q[i].x2 < x2min ? -1 : (q[i].x2 > x2max ? 1 : 0);
    ox3 = q[i].x3 < x3min ? -1 : (q[i].x3 > x3max ? 1 : 0);

    if (!pmy_block->pmy_mesh->multilevel) {
      fi1 = 0;
      fi2 = 0;
    } else {
      // reserved implementation for multilevel
    }

    int id = FindBufferID(ox1, ox2, ox3, fi1, fi2, pmy_block->pmy_mesh->maxneighbor_);

    if (id != -1) { // particle moved out of domain
      tmp = q[i];
      q[i] = q[j];
      q[j] = tmp;
      if (alive)
        bufid.push_back(id);
      else
        bufid.push_back(-1);
      j--;
    } else {  // particle in domain
      i++;
    }
  }
}
