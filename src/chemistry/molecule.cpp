// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "io_funcs.h"
#include "Molecule.h"
#include "Constants.h"
#include "utils.h"
#include "Parameter.h"

int Molecule::ntotal = 0;

std::ostream& operator<<(std::ostream &os, Molecule const& mol)
{
  os << "name: " << mol.m_name << std::endl
     << "molecular weight: " << mol.m_mu << " kg/mol" << std::endl
     << "heat capacity (cp): " << mol.m_cp << " J/(mol K)" << std::endl
     << "vaporization heat at triple point: " << mol.m_latent << " kJ/mol" << std::endl
     << "vapor pressure coefficients: " << mol.m_beta << " " << mol.m_gamma << std::endl
     << "standard entropy: " << mol.m_entropy << " J/(mol K)" << std::endl
     << "standard enthalpy: " << mol.m_enthalpy << " kJ/mol" << std::endl;
  return os;
}


Molecule::Molecule(std::string name = "") :
  name(name), mu(0), tr(0), pr(0), tc(0), pc(0), phase(GAS),
  m_cp(0), m_latent(0), m_entropy(0), m_enthalpy(0), m_gibbs(0),
  m_cliq(0), m_enliq(0), m_csld(0), m_ensld(0),
  m_beta(0), m_gamma(0), m_nshomate(0),
{
  prev = NULL;
  next = NULL;
  std::fill(m_shomate.begin(), m_shomate.end(), 0.);
  std::fill(m_shomate_sp.begin(), m_shomate_sp.end(), 0.);

  ntotal++;
}

Molecule* AddMolecule(std::string name)
{
  std::stringstream msg;
  Molecule *p = this;

  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule::AddMolecule. Molecule is empty. "
        << "Use new Molecle instred" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new Molecule(name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

void Molecule::Initialize(ParameterInput *pin)
{
  std::string gas = pin->GetString("chemistry", "gas");
  std::string cloud = pin->GetOrAddString("chemistry", "cloud", "");
  std::string folder = pin->GetString("chemistry", "folder");

  int ngas    = param.get_size("gas"),
      ncloud  = param.get_size("cloud");

    std::vector<std::string>
            agas(ngas),
            acloud(ncloud);

    std::string
            folder,
            name;

    param.get_entry("chem_folder", &folder);
    param.get_entry("gas", &agas.front(), &agas.back());
    param.get_entry("cloud", &acloud.front(), &acloud.back());

    for (int i = 0; i < ngas; ++i) {
        name = agas[i];
        auto molecule = std::make_shared<Molecule>(name);
        mlist.push_back(molecule);
        mlist.back()->load_chem_file(folder + name + ".chem");
        mlist.back()->set_to_simple();
    }

    for (int i = 0; i < ncloud; i++) {
        name = acloud[i];
        auto molecule = std::make_shared<Molecule>(name);
        mlist.push_back(molecule);
        std::string phase = name.substr(name.size() - 3, 3);
        auto gas = cfind_molecule(mlist, name.substr(0, name.size() - 3));
        // If cannot find an associated gas molecule, load a chemical file
        // without a phase name 
        if (gas == mlist.end()) {
            mlist.back()->load_chem_file(folder + name.substr(0, name.size() - 3) + ".chem");
            gas--;
        }
        if (phase == "(l)")
            mlist.back()->set_to_liquid(**gas);
        else
            mlist.back()->set_to_solid(**gas);
    }
}

void Molecule::LoadChemistryFile(std::string chemfile)
{
  std::stringstream inp(decomment(chemfile));
  Real junk;

  inp >> m_name >> m_mu
      >> m_entropy >> m_enthalpy >> m_gibbs >> m_cp
      >> m_tr >> m_pr
      >> m_tc >> m_pc; 
  inp >> m_nshomate;
  for (int i = 0; i < m_nshomate; i++) {
    inp >> m_shomate_sp.at(i) >> junk;
    for (int j = 0; j < NSHOMATE; j++)
      inp >> m_shomate.at(i * NSHOMATE + j);
  }
  m_shomate_sp.at(m_nshomate) = junk;

  inp >> m_cliq >> m_enliq
      >> m_csld >> m_ensld;

  m_mu *= 1.E-3;  // g/mol -> kg/mol
}

Real Molecule::Cp(Real T) const 
{
  int i = _locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  Real result;
  T /= 1.E3;
  result = m_shomate[i][0] + m_shomate[i][1]*T + m_shomate[i][2]*T*T + 
    m_shomate[i][3]*T*T*T + m_shomate[i][4]/(T*T);
  return result;
}

Real Molecule::Enthalpy(Real T) const
{
  int i = _locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  T /= 1.E3;
  return m_shomate[i][0]*T + 0.5 * m_shomate[i][1]*T*T +
    1./3. * m_shomate[i][2]*T*T*T + 1./4. * m_shomate[i][3]*T*T*T*T + 
    - m_shomate[i][4]/T + m_shomate[i][5];
}

Real Molecule::Entropy(Real T) const
{
  int i = _locate(m_shomate_sp.data(), T, m_nshomate + 1);
  if (i == m_nshomate) i = m_nshomate - 1; if (i < 0) i = 0;

  T /= 1.E3;
  return m_shomate[i][0] * log(T) + m_shomate[i][1]*T +
    0.5 * m_shomate[i][2]*T*T + 1./3. * m_shomate[i][3]*T*T*T + 
    - m_shomate[i][4]/(2.*T*T) + m_shomate[i][6];
}

Real Molecule::Latent(Real T) const {
  return m_latent + Enthalpy(T) - Enthalpy(m_tr) - m_cp * (T - m_tr) / 1.E3;
}

Real Molecule::SatVaporPres(Real T) const {
  return m_pr * exp((1. - m_tr / T) * m_beta - m_gamma * log(T / m_tr));
};

Real Molecule::SatVaporTemp(Real P, Real Tmin, Real Tmax, Real precision) const {
  Real Tsat, Tmin, Tmax;
  Tmin = 50.;
  Tmax = 2000.;
  int error = _root(Tmin, Tmax, precision, &Tsat, _SolveSatVaporTemp);
  std::stringstream msg;
  if (error) {
    msg << "SatVaporTemp failed" << std::endl
        << "Pressure = " << P << " Pa" << std::endl
        << "Tmin = " << Tmin << std::endl;
        << "Tmax = " << Tmax << std::endl;
        << "Saturation vapor pressure at Tmin = " << sat_vapor_p(Tmin) << std::endl;
        << "Saturation vapor pressure at Tmax = " << sat_vapor_p(Tmax) << std::endl;
        << *this << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return Tsat;
}

void Molecule::SetPhase(PhaseID pid) {
  if (pid == LIQUID) {
    m_name += "(l)";
    m_cp = m_cliq;
    m_latent = m_enliq;
    m_beta = (m_latent * 1.E3 - (cp(m_tr) - m_cp) * m_tr) / (Globals::Rgas * m_tr);
    m_gamma = (m_cp - cp(m_tr)) / Globals::Rgas;
    phase = LIQUID;
  } else if (pid == SOLID) {
    m_name = m_name + "(s)";
    m_cp = m_csld;
    m_latent = m_ensld;
    m_beta = (m_latent * 1.E3 - (cp(m_tr) - m_cp) * m_tr) / (Globals::Rgas * m_tr);
    m_gamma = (m_cp - cp(m_tr)) / Globals::Rgas;
  }
}

/*void Molecule::set_to_simple() {
    for (int i = 0; i < m_nshomate; i++)
        m_shomate.at(i * NSHOMATE) = m_cp;

    for (int i = 0; i < m_nshomate; i++)
        for (int j = 1; j < NSHOMATE; j++)
            m_shomate.at(i * NSHOMATE + j) = 0.;
}*/

Molecule* GetMolecule(std::string name)
{
  std::stringstream msg;
  Molecule *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule:GetMolecule " << name << " not found" <<
    std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return p;
}

Molecule const* GetMolecule(std::string name) const
{
  std::stringstream msg;
  Molecule const *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in Molecule:GetMolecule " << name << " not found" <<
    std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return p;
}

Real Molecule::_SolveSatVaporTemp(Real T, void *other)
{
  Real pres = *static_cast<Real*>(other);
  return SatVaporPres(T) - pres;
}
