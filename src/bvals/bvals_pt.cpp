//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_pt.cpp
//  \brief functions that apply BCs for particles

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "../particle/particle.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

void BoundaryValues::SendParticleBuffers(std::vector<Particle> const& pt, std::vector<int> const& bufid)
{
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;

  // clear send buffer
  for (int i = 0; i < pmb->pmy_mesh->maxneighbor_; ++i)
    particle_send_[i].clear();

  for (size_t i = 0; i < bufid.size(); ++i) {
    size_t j = pt.size() - i - 1;
    if (bufid[i] != -1) // if this particle is alive
      particle_send_[bufid[i]].push_back(pt[j]);
  }

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (nb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      particle_recv_[nb.targetid] = particle_send_[nb.bufid];
      pbl->pbval->particle_flag_[nb.targetid] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_ISend();
#endif
  }
}

bool BoundaryValues::ReceiveParticleBuffers(std::vector<Particle>& pt, std::vector<int>& bufid)
{
  MeshBlock *pmb = pmy_block_;
  bool flag = true;

  pt.resize(pt.size() - bufid.size());
  bufid.clear();

  // clear receiver buffer
  for (int i = 0; i < pmb->pmy_mesh->maxneighbor_; ++i)
    particle_recv_[i].clear();

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (particle_flag_[nb.bufid] == BNDRY_COMPLETED) continue;
    if (particle_flag_[nb.bufid] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe();
        MPI_IReceive();
        if (test == false) {
          flag = false;
          continue;
        }
        particle_flag_[nb.bufid] = BNDRY_ARRIVED;
      }
#endif
    }

    // the sign of nb.ox? might be flipped
    for (std::vector<Particle>::iterator it = particle_recv_[nb.bufid].begin();
        it != particle_recv_[nb.bufid].end(); ++it) {
      if (pmb->block_bcs[(nb.ox1+1)>>1] == PERIODIC_BNDRY) // 0:INNER_X1, 1:OUTER_X1
        it->x1 += nb.ox1*(pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min);

      if (pmb->block_bcs[2+(nb.ox2+1)>>1] == PERIODIC_BNDRY) // 2:INNER_X2, 3:OUTER_X2
        it->x2 += nb.ox2*(pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min);

      if (pmb->block_bcs[4+(nb.ox3+1)>>1] == PERIODIC_BNDRY) // 4:INNER_X3, 5:OUTER_X3
        it->x3 += nb.ox3*(pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min);

      // copy particle into the particle chain
      pt.push_back(*it);
    }
  }

  return flag;
}
