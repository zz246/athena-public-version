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

  for (size_t i = 0; i < bufid.size(); ++i) {
    size_t j = pt.size() - i - 1;
    if (bufid[i] != -1) // if this particle is alive
      particle_send_[bufid[i]].push_back(pt[j]);
  }

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (nb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(nb.gid);
      hydro_recv_[nb.targetid] = hydro_send_[nb.bufid];
      pbl->pbval->particle_flag_[nb.targetid] = BNDY_ARRIVED;
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

  for (int n = 0; n < pmb->nneighbor; ++n) {
    NeighborBlock &nb = pmb->neighbor[n];

    if (particle_flag_[nb.bufid] == BNDRY_COMPLETED) continue;
    if (particle_flag_[nb.bufid] == BNDRY_WAITING) {
      if (nb.rank == Globals::my_rank) { // on the same process
        flag = false;
        continue
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

    // remember special considerations for periodic boundary condition

    for (size_t i = 0; i < hydro_recv_[nb.bufid].size(); ++i)
      pt.push_back(hydro_recv_[i]);
    hydro_recv_[nb.bufid].clear();
  }

  return flag;
}
