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

void BoundaryValues::DetachParticle(
  std::vector<Particle> const& pt, std::vector<int> const& bufid, int pid)
{
  for (size_t i = 0; i < bufid.size(); ++i)
    if (bufid[i] >= 0) { // if particle is still alive
      size_t j = pt.size() - i - 1;
      particle_send_[bufid[i]].push_back(pt[j]);
    }
  for (int i = 0; i < pmy_block_->pmy_mesh->maxneighbor_; ++i)
    particle_num_send_[i][pid] = particle_send_[i].size();
}

void BoundaryValues::SendParticleBuffers()
{
  MeshBlock *pmb = pmy_block_;

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (nb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      pbl->pbval->particle_recv_[nb.targetid] = particle_send_[nb.bufid];
      for (int i = 0; i < ParticleGroup::ntotal; ++i)
        pbl->pbval->particle_num_recv_[nb.targetid][i] = particle_num_send_[nb.bufid][i];
      pbl->pbval->particle_flag_[nb.targetid] = BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else { // MPI
      MPI_Start(&req_particle_num_send_[nb.bufid]);
      int tag = CreateBvalsMPITag(nb.lid, TAG_PARTICLE_NUM, nb.targetid);
      int ssize = particle_send_[nb.bufid].size();
      MPI_Isend(particle_send_[nb.bufid].data(), ssize, MPI_PARTICLE,
                nb.rank, tag, MPI_COMM_WORLD, &req_particle_send_[nb.bufid]);
    }
#endif
  }
}

bool BoundaryValues::ReceiveParticleBuffers(
  std::vector<Particle>& pt, std::vector<int>& bufid, int pid)
{
  MeshBlock *pmb = pmy_block_;
  bool flag = true;
  int test, rsize, tag;

  // clear 
  pt.resize(pt.size() - bufid.size());
  bufid.clear();

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (particle_flag_[nb.bufid] == BNDRY_COMPLETED) continue;
    if (particle_flag_[nb.bufid] == BNDRY_INIT) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        MPI_Test(&req_particle_num_recv_[nb.bufid], &test, MPI_STATUS_IGNORE);
        if (test == false) {
          flag = false;
          continue;
        } else {
          rsize = particle_num_recv_[nb.bufid][ParticleGroup::ntotal - 1];
          particle_recv_[nb.bufid].resize(rsize);
          tag = CreateBvalsMPITag(nb.lid, TAG_PARTICLE, nb.bufid);
          MPI_Irecv(particle_recv_[nb.bufid].data(), rsize, MPI_PARTICLE,
                    nb.rank, tag, MPI_COMM_WORLD, &req_particle_recv_[nb.bufid]);
          particle_flag_[nb.bufid] = BNDRY_WAITING;
        }
      }
#endif
    }
    if (particle_flag_[nb.bufid] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        MPI_Test(&req_particle_recv_[nb.bufid], &test, MPI_STATUS_IGNORE);
        if (test == false) {
          flag = false;
          continue;
        }
        particle_flag_[nb.bufid] = BNDRY_ARRIVED;
      }
#endif
    }
    // attach particle to the particle chain
    std::vector<Particle>::iterator begin = particle_recv_[nb.bufid].begin();
    std::vector<Particle>::iterator it = pid == 0 ? begin : begin + particle_num_recv_[nb.bufid][pid - 1];
    for (; it != begin + particle_num_recv_[nb.bufid][pid]; ++it) {
      if (pmb->block_bcs[(nb.ox1+1)>>1] == PERIODIC_BNDRY) // 0:INNER_X1, 1:OUTER_X1
        it->x1 += nb.ox1*(pmb->pmy_mesh->mesh_size.x1max - pmb->pmy_mesh->mesh_size.x1min);

      if (pmb->block_bcs[2+(nb.ox2+1)>>1] == PERIODIC_BNDRY) // 2:INNER_X2, 3:OUTER_X2
        it->x2 += nb.ox2*(pmb->pmy_mesh->mesh_size.x2max - pmb->pmy_mesh->mesh_size.x2min);

      if (pmb->block_bcs[4+(nb.ox3+1)>>1] == PERIODIC_BNDRY) // 4:INNER_X3, 5:OUTER_X3
        it->x3 += nb.ox3*(pmb->pmy_mesh->mesh_size.x3max - pmb->pmy_mesh->mesh_size.x3min);

      if (it->x1 < pmb->block_size.x1min || it->x1 > pmb->block_size.x1max
          || it->x2 < pmb->block_size.x2min || it->x2 > pmb->block_size.x2max
          || it->x3 < pmb->block_size.x3min | it->x3 > pmb->block_size.x3max) {
        flag = false;
        std::stringstream msg;
        msg << "### FATAL ERROR in ReceiveParticleBuffers: particle moved out of MeshBlock limits" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      pt.push_back(*it);
    }

    if (pid == ParticleGroup::ntotal - 1)
      particle_flag_[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  return flag;
}
