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

  MeshBlock *pmb = pmy_block;
  Hydro *phydro = pmb->phydro;

  v1.InitWithShallowSlice(phydro->w,4,IM1,1);
  v2.InitWithShallowSlice(phydro->w,4,IM2,1);
  v3.InitWithShallowSlice(phydro->w,4,IM3,1);

  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

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

    bool alive = particle_fn_(pmb, q[i], time, dt);

    // take care of reflective boundary condition
    if (q[i].x1 < x1min && pmb->block_bcs[0] == REFLECTING_BNDRY)
      q[i].x1 = 2*x1min - q[i].x1;
    if (q[i].x1 > x1max && pmb->block_bcs[1] == REFLECTING_BNDRY)
      q[i].x1 = 2*x1max - q[i].x1;
    if (q[i].x2 < x2min && pmb->block_bcs[2] == REFLECTING_BNDRY)
      q[i].x2 = 2*x2min - q[i].x2;
    if (q[i].x2 > x2max && pmb->block_bcs[3] == REFLECTING_BNDRY)
      q[i].x2 = 2*x2max - q[i].x2;
    if (q[i].x3 < x3min && pmb->block_bcs[4] == REFLECTING_BNDRY)
      q[i].x3 = 2*x3min - q[i].x3;
    if (q[i].x3 > x3max && pmb->block_bcs[5] == REFLECTING_BNDRY)
      q[i].x3 = 2*x3max - q[i].x3;

    ox1 = q[i].x1 < x1min ? -1 : (q[i].x1 > x1max ? 1 : 0);
    ox2 = q[i].x2 < x2min ? -1 : (q[i].x2 > x2max ? 1 : 0);
    ox3 = q[i].x3 < x3min ? -1 : (q[i].x3 > x3max ? 1 : 0);

    if (!pmb->pmy_mesh->multilevel) {
      fi1 = 0;
      fi2 = 0;
    } else {
      // reserved implementation for multilevel
    }

    int id = FindBufferID(ox1, ox2, ox3, fi1, fi2, pmb->pmy_mesh->maxneighbor_);

    if (alive && (id == -1)) { // particle is alive and inside domain
      i++;
    } else {  // particle deseased or moved out of the domain
      tmp = q[i];
      q[i] = q[j - 1];
      q[j - 1] = tmp;
      bufid.push_back(alive ? id : -1); // Note that bufid is reversed
      j--;
    }
  }
}
