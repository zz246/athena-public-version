//! \file sw.cpp
//  \brief global shallow water model

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
}
