//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      File name:     G4NRFNuclearLevel
//      Author:        David Jordan (david.jordan@pnl.gov)
//      Creation date: October 2006
//
//      Adapted from:  G4NuclearLevel,
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//      Creation date: 24 October 1998
//
//      Modifications:
//        14 Jan. 2015, Jayson Vavrek (jvavrek@mit.edu)
//              Moderate to major bug fixes plus enhanced documentation.
//
//        14 Oct. 2006, David Jordan (DavidJordan@pnl.gov)
//              Adapted G4NuclearLevel class to NRF physics process.
//
//        G4NuclearLevel history follows:
//
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added IC probability when calculate the channel probabilities in
//              MakeProbabilities().
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions.
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//
// -------------------------------------------------------------------

// -------------------------------------------------------------------
//
// Description:
//
// G4NRFNuclearLevel is the class used to instantiate a particular level of a
// particular nucleus; its members include the energy, spin, etc of the level.
// To build the levels, data is imported using the G4NRFNuclearLevelManager
// class and an attempt to coordinate the different databases is made here,
// since various databases only compile values for certain quantities (e.g.)
// the energy levels and spins but not multipole moments. To match up levels
// across databases, the level energies are used as identifiers; invariably,
// the energies differ slightly, so a tolerance of 10.0 eV is specified by
// default.
//
//
// Updates by Jayson Vavrek:
//
// 0) Modernized (2015) C++ standards: #include's, using's; also
//    G4SystemOfUnits.hh and G4PhysicalConstants.hh
//
// 1) In RefreshGammas(), implemented a workaround for erroneous large
//    negative energy values in databases
//
// 2) In RefreshGammas(), fixed variable declaration placement when looking
//    for Ediff_min
//
// 3) In RefreshWidth(), implemented the same workaround as in RefreshGammas()
//
// 4) In RefreshWidth(), implemented a fix for levels with empirically-unknown
//    spins. If the spin of the level is unknown, it would previously get set to
//    -999 and cause an exit() error. Now, by using Get/SetInvalidLevel() of
//    G4NRFNuclearLevel, we can mark the level as invalid, then in
//    G4NRFNuclearLevelManager, we use the function delete_bad_spins() to erase
//    such levels in their entirety.
//
// 5) Addition of standalone NRF gamma database generation. See PrintAllTabular()
//    or the G4NRF.cc source file.
//
// -------------------------------------------------------------------

#include "G4NRFNuclearLevel.hh"

#include <fstream>

#include "iomanip"
using std::setw;

#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4NRF.hh"


const G4double EDIFF_TOL_gamma = 10.0*eV; // default is 10.0*eV, probably good enough
const G4double EDIFF_TOL_level = 10.0*eV;

G4int G4NRFNuclearLevel::Increment(G4int aF) {
  static G4int instanceCount = 0;
  instanceCount += aF;
  return instanceCount;
}

G4NRFNuclearLevel::G4NRFNuclearLevel(const G4int nLevel,
             const G4int nucleusZ, const G4int nucleusA,
             const G4double E_1stExcitedState,
             const G4double energy, const G4double halfLife,
             const G4double angularMomentum,
             const std::vector<double>& eGamma,
             const std::vector<double>& wGamma,
             const std::vector<double>& polarities,
             const std::vector<double>& kCC, const std::vector<double>& l1CC,
             const std::vector<double>& l2CC, const std::vector<double>& l3CC,
             const std::vector<double>& m1CC, const std::vector<double>& m2CC,
             const std::vector<double>& m3CC, const std::vector<double>& m4CC,
             const std::vector<double>& m5CC, const std::vector<double>& nPlusCC,
             const std::vector<double>& totalCC,
             const G4bool Verbose
             ) {
  _nLevel            = nLevel;
  _nucleusZ          = nucleusZ;
  _nucleusA          = nucleusA;
  _E_1stExcitedState = E_1stExcitedState;
  _energy            = energy;
  _halfLife          = halfLife;
  _angularMomentum   = angularMomentum;
  _Verbose = Verbose;

  _cross_sec_interp_func = NULL;

  invalidLevel = false;

  unsigned int i;
  for (i = 0; i < eGamma.size(); i++) {
    _energies.push_back(eGamma[i]);
    _weights.push_back(wGamma[i]);
    _polarities.push_back(polarities[i]);
    _kCC.push_back(kCC[i]);
    _l1CC.push_back(l1CC[i]);
    _l2CC.push_back(l2CC[i]);
    _l3CC.push_back(l3CC[i]);
    _m1CC.push_back(m1CC[i]);
    _m2CC.push_back(m2CC[i]);
    _m3CC.push_back(m3CC[i]);
    _m4CC.push_back(m4CC[i]);
    _m5CC.push_back(m5CC[i]);
    _nPlusCC.push_back(nPlusCC[i]);
    _totalCC.push_back(totalCC[i]);
  }
  _nGammas = _energies.size();

  //  G4cout << "_nGammas: " << _nGammas << G4endl;

  MakeProbabilities();
  MakeCumProb();

  _Tau = hbar_Planck * log(2)/_halfLife;

  // Calls to member functions for NRF-specific information follow.
  // "Refresh" refers to augmenting level information read from
  // default G4 photon evaporation database.

  RefreshWidth();   // update level width information
  RefreshGammas();  // update gamma multipolarity information
  MakeWidth0();     // set up g.s. transition partial width

//#endif

  if (_Verbose) {
    PrintAll();
  }
}


G4NRFNuclearLevel::~G4NRFNuclearLevel()
{;}


G4bool G4NRFNuclearLevel::operator==(const G4NRFNuclearLevel &right) const {
  return (this == (G4NRFNuclearLevel *) &right);
}


G4bool G4NRFNuclearLevel::operator!=(const G4NRFNuclearLevel &right) const {
  return (this != (G4NRFNuclearLevel *) &right);
}


G4bool G4NRFNuclearLevel::operator<(const G4NRFNuclearLevel &right) const {
  if (_energy < right.Energy()) {
    return true;
  } else {
    return false;
  }
}


const std::vector<double>& G4NRFNuclearLevel::GammaEnergies() const {
  return _energies;
}

const std::vector<double>& G4NRFNuclearLevel::GammaWeights() const {
  return _weights;
}

const std::vector<double>& G4NRFNuclearLevel::GammaProbabilities() const {
  return _prob;
}

const std::vector<double>& G4NRFNuclearLevel::GammaCumulativeProbabilities() const {
  return _cumProb;
}

const std::vector<double>& G4NRFNuclearLevel::GammaPolarities() const {
  return _polarities;
}

const std::vector<double>& G4NRFNuclearLevel::KConvertionProbabilities() const {
  return _kCC;
}

const std::vector<double>& G4NRFNuclearLevel::L1ConvertionProbabilities() const {
  return _l1CC;
}

const std::vector<double>& G4NRFNuclearLevel::L2ConvertionProbabilities() const {
  return _l2CC;
}

const std::vector<double>& G4NRFNuclearLevel::L3ConvertionProbabilities() const {
  return _l3CC;
}

const std::vector<double>& G4NRFNuclearLevel::M1ConvertionProbabilities() const {
  return _m1CC;
}

const std::vector<double>& G4NRFNuclearLevel::M2ConvertionProbabilities() const {
  return _m2CC;
}

const std::vector<double>& G4NRFNuclearLevel::M3ConvertionProbabilities() const {
  return _m3CC;
}

const std::vector<double>& G4NRFNuclearLevel::M4ConvertionProbabilities() const {
  return _m4CC;
}

const std::vector<double>& G4NRFNuclearLevel::M5ConvertionProbabilities() const {
  return _m5CC;
}

const std::vector<double>& G4NRFNuclearLevel::NPlusConvertionProbabilities() const {
  return _nPlusCC;
}

const std::vector<double>& G4NRFNuclearLevel::TotalConvertionProbabilities() const {
  return _totalCC;
}

const std::vector<int>& G4NRFNuclearLevel::MultipoleNumModes() const {
  return _Num_multipole;
}

const std::vector<char>& G4NRFNuclearLevel::MultipoleMode1() const {
  return _Multipole_mode1;
}

const std::vector<int>& G4NRFNuclearLevel::MultipoleL1() const {
  return _Multipole_L1;
}

const std::vector<char>& G4NRFNuclearLevel::MultipoleMode2() const {
  return _Multipole_mode2;
}

const std::vector<int>& G4NRFNuclearLevel::MultipoleL2() const {
  return _Multipole_L2;
}

const std::vector<double>& G4NRFNuclearLevel::MultipoleMixingRatio() const {
  return _Multipole_mixing_ratio;
}

const std::vector<int>& G4NRFNuclearLevel::MultipoleMixRatioSignFlag() const {
  return _Multipole_mixing_sign_flag;
}

G4int G4NRFNuclearLevel::nLevel() const {
  return _nLevel;
}

G4int G4NRFNuclearLevel::Z() const {
  return _nucleusZ;
}

G4int G4NRFNuclearLevel::A() const {
  return _nucleusA;
}

G4double G4NRFNuclearLevel::AngularMomentum() const {
  return _angularMomentum;
}

G4double G4NRFNuclearLevel::Parity() const {
  return _parity;
}

G4double G4NRFNuclearLevel::HalfLife() const {
  return _halfLife;
}

G4int G4NRFNuclearLevel::NumberOfGammas() const {
  return _nGammas;
}

G4double G4NRFNuclearLevel::Width() const {
  return _Tau;
}

G4double G4NRFNuclearLevel::Width0() const {
  return _Tau0;
}

G4double G4NRFNuclearLevel::MaxGammaEnergy() const {
  if (_nGammas > 0)
    return _energies[_nGammas-1];
  else
    return 0.0;
}

G4bool G4NRFNuclearLevel::GetInvalidLevel() {
  return invalidLevel;
}

void G4NRFNuclearLevel::SetInvalidLevel() {
  invalidLevel = true;
}

G4double G4NRFNuclearLevel::SelectGamma(G4int& igamma) const {
  // Returns gammaEnergy = 0 if no gammas
  // gammaEnergy > 0.0 if gamma transition
  // gammaEnergy < 0.0 if conversion electron
  // Pass-by-reference parameter igamma is the
  // index of gamma selected [0, 1, ...,
  // _nGammas-1].

  G4double gammaEnergy = 0.0;

  if (_nGammas > 0) {
    if ((_nGammas == 1) && (_weights[0] == 0.0)) {
      gammaEnergy = -_energies[0];
      igamma = 0;
    } else {
      G4double random = G4UniformRand();

      G4int iGamma = 0;
      for (iGamma = 0; iGamma < _nGammas; iGamma++) {
        if (random <= _cumProb[iGamma]) break;
      }

      gammaEnergy = _energies[iGamma];
      igamma = iGamma;

      // now decide whether Internal Coversion electron should be emitted instead
      random = G4UniformRand();
      G4double tot_cc = _totalCC[iGamma];
      G4double IC_prob = tot_cc/(tot_cc + 1.0);

      if (random <= IC_prob)
        gammaEnergy *= -1.0;
    }
  }

  return gammaEnergy;
}


void G4NRFNuclearLevel::PrintAll() const {
  G4cout << "---- Level energy (MeV): " << _energy/MeV << ", angular momentum: "
         << _angularMomentum << ", half life (ps): " << _halfLife/picosecond
         << ", #photons: " << _nGammas << G4endl
         << "width (eV): " << _Tau/eV
         << ", width0 (eV): " << _Tau0/eV
         << G4endl;
  G4cout << "    E_1stExcitedState (MeV): " << _E_1stExcitedState/MeV << G4endl;

  G4int i;
  G4cout << "     Gammas: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _energies[i] << " "; }

  G4cout << G4endl << "     Weights: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _weights[i] << " "; }

  G4cout << G4endl << "     Relative transition probabilities ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _prob[i] << " "; }

  G4cout << G4endl << "     Cumulative probabilities: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _cumProb[i] << " "; }

  G4cout << G4endl << "     Polarities: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _polarities[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     NumMultipoles: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Num_multipole[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleMode1: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_mode1[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleL1: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_L1[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleMode2: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_mode2[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleL2: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_L2[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleMixingRatio: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_mixing_ratio[i] << " "; }

  G4cout << G4endl;
  G4cout << G4endl << "     MultipoleMixRatio SignFlag: ";
  for (i = 0; i < _nGammas; ++i) { G4cout << _Multipole_mixing_sign_flag[i] << " "; }
  G4cout << G4endl;

  return;
}

// print relevant info for standalone code in tabular format
// hacky way to pass gsSpin and TDebye values known to NuclearLevelManager but not to NuclearLevel
void G4NRFNuclearLevel::PrintAllTabular(ofstream& file, double gsSpin, double TDebye) const {
  for (int i = 0; i < _nGammas; i++) {
    file << _nucleusZ << " "
         << _nucleusA << " "
         << _energy/MeV << " "
         << _energies[i]/MeV << " "
         << _Tau/MeV << " "
         << _prob[i] << " "
         << GetGSProb() << " "
         << gsSpin << " "
         << _angularMomentum << " "
         << TDebye << " "
         << G4endl;
  }
}


G4double G4NRFNuclearLevel::GetGSProb() const {
  double GSprob = 0.0;
  for (int i = 0; i < _nGammas; ++i) {
    if (fabs(_energies[i] - _energy) <= EDIFF_TOL_level) {
      GSprob = _prob[i];
    }
  }
  return GSprob;
}


void G4NRFNuclearLevel::MakeProbabilities() {
  G4double sum = 0.;
  G4int i = 0;
  for (i = 0; i < _nGammas; ++i) {
    sum += _weights[i]*(1+_totalCC[i]);
  }

  for (i = 0; i < _nGammas; ++i) {
    if (sum > 0.) {
      _prob.push_back(_weights[i]*(1+_totalCC[i])/ sum);
    } else {
      _prob.push_back(1./_nGammas);
    }
  }
  return;
}


void G4NRFNuclearLevel::MakeCumProb() {
  if (_nGammas > 0) {
    G4double sum = _prob[0];
    _cumProb.push_back(sum);

    G4int i = 0;
    for (i = 1; i < _nGammas; i++) {
      sum += _prob[i];
      _cumProb.push_back(sum);
    }
  }
  return;
}

void G4NRFNuclearLevel::RefreshGammas() {
  // NRF-specific routine
  //
  // Supplement gamma info with multipolarity information.  Uses the tables
  //        gamma_table_nnn.dat
  // in the directory pointed to by $G4NRFGAMMADATA.
  //
  // These tables have the format
  //    Z   E_level(keV)   Num_Gammas
  //       E_gamma1(keV)  Num_Multipoles  [M/E  L1   M/E   L2   delta   sign_flag]
  //       E_gamma2(keV)  Num_Multipoles  [M/E  L1   M/E   L2   delta   sign_flag]
  //       E_gamma3(keV)  Num_Multipoles  [M/E  L1   M/E   L2   delta   sign_flag]
  //       ....
  //       E_gammaN(keV)  Num_Multipoles  [M/E  L1   M/E   L2   delta   sign_flag]
  //
  // where there are i = 1, 2, ..., N = Num_Gammas gamma records in each (Z,E_level)
  // block.  The number of multipoles for a particular gamma can be 0, 1, or 2.  If
  // num_multipoles = 1, then the multipole mode and angular momentum are read from
  // the next two entries in the line.  The mode is specified by the single character
  // "E", "M", "D" (dipole), or "Q" (quadrupole).  The angular momentum L is an
  // integer.  If num_multipoles = 2, a second set of (mode, L) values is read, in
  // addition to the mixing ratio, delta, and a sign flag (indicates whether the
  // sign of the mixing ratio is known.)


  // Open gamma table file

  char* env = getenv("G4NRFGAMMADATA");
  if (!env) {
    G4cout << "G4NRFNuclearLevel: please set the G4NRFGAMMADATA environment variable" << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(31);
  }

  G4String dir = G4String(env) + '/';

  std::ostringstream streamName;

  streamName << "gamma_table_";
  if (_nucleusA < 10) {
    streamName << "00" << _nucleusA;
  } else if (_nucleusA < 100) {
    streamName << "0"  << _nucleusA;
  } else {
    streamName << _nucleusA;
  }

  streamName << ".dat";

  G4String gamma_table_filename = dir + G4String(streamName.str());

  if (_Verbose) {
    G4cout << "**************************** In G4NRFNuclearLevel::RefreshGammas()." << G4endl;
    G4cout << "E (keV): " << _energy/keV << G4endl
     << " Z: " << _nucleusZ << G4endl
     << " gamma_table_file: " << gamma_table_filename << G4endl;
  }

  std::ifstream GammaTableFile(gamma_table_filename, std::ios::in);
  if (!GammaTableFile) {
    G4cout << "Error in G4NRFNuclearLevel::RefreshGammas()" << G4endl;
    G4cout << "Could not open file "<< gamma_table_filename << G4endl;
    G4cout << "Expected to find this file in directory" << G4endl;
    G4cout << dir << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(32);
  }


  // Locate level in the gamma table
  G4bool found_level = false;
  char line[130];

  G4int Z;
  G4int num_gammas;
  G4double E_level;
  G4double closest_level  = 0.0*keV;
  G4double min_level_diff = 1000.0*keV;

  // level-matching algorithm, part 1
  while ((!found_level) && GammaTableFile.getline(line, 130)) {
    std::istringstream line_str(line);
    line_str >> Z >> E_level >> num_gammas;

    // deal with the bizarre lines with +X in NNDC
    if (E_level <= -2.4e6) {
      G4double E_init = E_level;
      E_level *= -1.0;
      E_level -= 2.4e6;
    }

    if (Z == _nucleusZ) {
      E_level *= keV;
      G4double Ediff = E_level - _energy;
      if (Ediff < 0.0) Ediff *= -1.0;
      if (Ediff < min_level_diff) {
        min_level_diff = Ediff;
        closest_level = E_level;
      }
      if (Ediff <= EDIFF_TOL_level) {
        found_level = true; // only place found_level is set in RefreshGammas();
        if (_Verbose) G4cout << "RefreshGammas: Found Z = " << Z << " level at E= " << E_level/keV << G4endl;
      }
    }

    if (!found_level) {
      for (G4int iskip = 0; iskip < num_gammas; iskip++) {
        GammaTableFile.getline(line, 130);
      }
    }
  } // end of level-matching alogrithm while loop


  // the algorithm looks for an E_level in the file to match the required _energy
  // if one is found, found_level = true and the algorithm proceeds past this next check
  // if not, abort by exit(1)

  if (!found_level) {
    if (false) { // may be useful later
      G4cout << "Error in G4NRFNuclearLevel::RefreshGammas." << G4endl;
      G4cout << "RefreshGammas: Could not find level " << _energy/keV << " keV." << G4endl;
      G4cout << "Closest level " << closest_level/keV << " keV, off by " << min_level_diff/keV << " keV." << G4endl;
      G4cout << "Nucleus A: " << _nucleusA << " Z: " << _nucleusZ << G4endl;
      G4cout << "Table file: " << gamma_table_filename << G4endl;
      G4cout << "Setting invalid level at E/MeV = " << _energy/MeV << G4endl;
      G4cout << G4endl;
    }
    this->SetInvalidLevel();
  }

  // Store the list of (gamma, multipole info) for the level
  const    G4int MAX_GAMMAS = 1000;
  G4double Egamma_arr[MAX_GAMMAS];
  G4int    num_multipole_arr[MAX_GAMMAS];
  char     multipole_mode1_arr[MAX_GAMMAS];
  G4int    multipole_L1_arr[MAX_GAMMAS];
  char     multipole_mode2_arr[MAX_GAMMAS];
  G4int    multipole_L2_arr[MAX_GAMMAS];
  G4double mixing_ratio_arr[MAX_GAMMAS];
  G4int    mixing_ratio_sign_flag_arr[MAX_GAMMAS];

  const    G4int MAX_MULT = 10;
  char     multipole_mode_arr[MAX_MULT];
  G4int    multipole_L_arr[MAX_MULT];

  G4int    jgamma;

  for (G4int jgamma = 0; jgamma < num_gammas; jgamma++) {
    GammaTableFile >> Egamma_arr[jgamma] >> num_multipole_arr[jgamma];
    Egamma_arr[jgamma] *= keV;
    if (num_multipole_arr[jgamma] > 0) {
      for (G4int imult = 0; imult < num_multipole_arr[jgamma]; imult++) {
        GammaTableFile >> multipole_mode_arr[imult] >> multipole_L_arr[imult];
      }
      multipole_mode1_arr[jgamma] = multipole_mode_arr[0];
      multipole_L1_arr[jgamma]    = multipole_L_arr[0];
      if (num_multipole_arr[jgamma] > 1) {
        multipole_mode2_arr[jgamma] = multipole_mode_arr[1];
        multipole_L2_arr[jgamma]    = multipole_L_arr[1];
        GammaTableFile >> mixing_ratio_arr[jgamma] >> mixing_ratio_sign_flag_arr[jgamma];
      }
    }
  }


  // Loop over the level's gammas and attempt to match them
  // up with the gammas stored from the multipole-info table.
  // If found, the multipole info for the gamma is stored
  // in the corresponding member variable vector.  If not
  // found, the level is marked as invalid.
  G4int i;

  if (_Verbose) {
    G4cout << "G4NRFNuclearLevel::RefreshGammas -- " << G4endl;
    G4cout << "Beginning update loop for mixing ratio info." << G4endl;
    G4cout << "_nGammas: " << _nGammas << G4endl;
  }

  // level-matching algorithm, part 2
  for (i = 0; i < _nGammas; i++) {
    G4double E_g = _energies[i];

    G4double Ediff_min = DBL_MAX; // fixed the placements of these (were outside of the for loop before)
    G4int jgamma_min = 0;

    if (_Verbose)
      G4cout << "  Trying to locate gamma of energy E_g (keV) = " << E_g/keV << G4endl;

    G4bool found_gamma = false;
    jgamma = 0;
    while ((!found_gamma) && (jgamma < num_gammas)) {
      G4double Ediff = fabs(E_g - Egamma_arr[jgamma]);
      if (Ediff < Ediff_min) {
        Ediff_min = Ediff;
        jgamma_min = jgamma;
      }

      // added relative tolerance of 0.1% - solves gamma errors in Pu but not level errors
      if (Ediff <= EDIFF_TOL_gamma || Ediff/E_g <= 1.0e-3) {
        found_gamma = true; // only place found_gamma is set in Refresh_Gammas()

        if (_Verbose) G4cout << "RefreshGammas: Found gamma (keV): " << E_g/keV << G4endl;
      } else {
        jgamma++;
      }
    } // end while loop


    // All's well, located gamma in multipole-properties table.
    // Store multipole info for this gamma, distinguishing
    // cases where there are 0, 1, and 2 different multipole
    // modes.  Note: To keep push_back ordering in step with
    // gamma vector, make sure to store dummy values if the
    // info is not available/relevant for the current gamma.
    if (found_gamma) {
      G4int nmult = num_multipole_arr[jgamma];

      _Num_multipole.push_back(nmult);

      if (nmult > 0) {
        char mode1 = multipole_mode1_arr[jgamma];
        G4int L1 = multipole_L1_arr[jgamma];
        _Multipole_mode1.push_back(mode1);
        _Multipole_L1.push_back(L1);

        if (nmult > 1) {
          char mode2 = multipole_mode2_arr[jgamma];
          G4int L2 = multipole_L2_arr[jgamma];
          G4double mixing_ratio = mixing_ratio_arr[jgamma];
          G4int sign_flag = mixing_ratio_sign_flag_arr[jgamma];
          _Multipole_mode2.push_back(mode2);
          _Multipole_L2.push_back(L2);
          _Multipole_mixing_ratio.push_back(mixing_ratio);
          _Multipole_mixing_sign_flag.push_back(sign_flag);
        } else {
          _Multipole_mode2.push_back('X');
          _Multipole_L2.push_back(-1);
          _Multipole_mixing_ratio.push_back(-999.0);
          _Multipole_mixing_sign_flag.push_back(-1);
        }
      } else {
        _Multipole_mode1.push_back('X');
        _Multipole_L1.push_back(-1);
        _Multipole_mode2.push_back('X');
        _Multipole_L2.push_back(-1);
        _Multipole_mixing_ratio.push_back(-999.0);
        _Multipole_mixing_sign_flag.push_back(-1);
      }
    } else {
     if (false) { // may be useful later
      G4cout << "Error in G4NRFNuclearLevel::RefreshGammas." << G4endl;
      G4cout << "Could not locate gamma in table." << G4endl;
      G4cout << "Looking for gamma (keV): " << E_g/keV << G4endl;
      G4cout << "Nucleus A: " << _nucleusA << " Z: " << _nucleusZ << G4endl;
      G4cout << "From level E (keV) = " << _energy/keV << G4endl;

      G4cout << "There are " << num_gammas << " candidates in table "
             << gamma_table_filename << " as follows: " << G4endl;
      for (jgamma = 0; jgamma < num_gammas; jgamma++)
        G4cout << jgamma << " " << Egamma_arr[jgamma]/keV << G4endl;

      G4cout << "Minimum difference occurs for gamma #" << jgamma_min << G4endl;
      G4cout << "Ediff_min (keV): " << Ediff_min/keV << G4endl;
      //G4cout << "Setting invalid level at E = " << _energy << G4endl;
      //G4cout << G4endl;
     }
     //G4cout << "Setting invalid level at E = " << _energy << G4endl;
     this->SetInvalidLevel();
    } // found gamma in the multipole table?
  } // loop over gammas in this level
}

void G4NRFNuclearLevel::RefreshWidth() {
  // NRF-specific routine
  //
  // Supplement gamma info from photon cascade database
  // with level width information.  Draws upon files
  //        level_table_nnn.dat
  // in the directory pointed to by $G4NRFGAMMADATA.


  char* env = getenv("G4NRFGAMMADATA");
  if (!env) {
    G4cout << "G4NRFNuclearLevel: please set the G4NRFGAMMADATA environment variable" << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(35);
  }

  G4String dir = G4String(env) + '/';

  std::ostringstream streamName;

  streamName << "level_table_";
  if (_nucleusA < 10) {
    streamName << "00" << _nucleusA;
  } else if (_nucleusA < 100) {
    streamName << "0"  << _nucleusA;
  } else {
    streamName << _nucleusA;
  }

  streamName << ".dat";

  G4String levels_filename = dir + G4String(streamName.str());

  if (_Verbose) {
    G4cout << "********************************************** In G4NRFNuclearLevel::RefreshWidth()." << G4endl;
    G4cout << "E (keV): " << _energy/keV << " levels_file: " << levels_filename << G4endl;
  }

  std::ifstream LevelsFile(levels_filename, std::ios::in);
  if (!LevelsFile) {
    G4cout << "Error in G4NRFNuclearLevel::RefreshWidth()" << G4endl;
    G4cout << "Could not open file "<< levels_filename << G4endl;
    G4cout << "Expected to find this file in directory" << G4endl;
    G4cout << dir << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(36);
  }

  G4bool found_level = false;

  char line[130];

  G4int Z, A;
  G4double E_level;
  G4double spin, parity;
  G4double T_half, E_width;

  // level-matching algorithm, part 3
  while ((!found_level) && (!LevelsFile.eof())) {
    LevelsFile.getline(line, 130);
    std::istringstream line_str(line);
    line_str >> Z >> A >> E_level >> spin >> parity
       >>  T_half >> E_width >> _Ewidth_gamma
       >> _Ewidth_gamma0 >> _Ewidth_p
       >> _Ewidth_n >> _Ewidth_alpha;

    // deal with the bizarre lines with +X in NNDC
    if (E_level <= -2.4e6) {
      G4double E_init = E_level;
      E_level *= -1.0;
      E_level -= 2.4e6;
    }

    if ((Z == _nucleusZ) && (A == _nucleusA)) {
      E_level *= keV;
      G4double Ediff = E_level - _energy;
      if (Ediff < 0.0) Ediff *= -1.0;

      // the actual check on energy matching; we can be less strict if we know spins match
      if (Ediff <= EDIFF_TOL_level) {
        //G4cout << "spin = " << spin << " am = " << _angularMomentum << G4endl;
        found_level = true; // only place found_level is set in RefreshWidth()
        E_width *= eV;
        T_half  *= second;

        if (_Verbose) {
          G4cout << "Replacing _Tau (eV) = "
                 << _Tau/eV
                 << " with Ewidth (eV) = "
                 << E_width/eV << G4endl;
          G4cout << "Replacing _halfLife (sec) = "
                 << _halfLife/second
                 << " with T_half (sec) = "
                 << T_half/second << G4endl;
        }

        // deal with levels that have empirically-unknown spins
        if (spin == -999) this->SetInvalidLevel();

        if (spin != _angularMomentum && !this->GetInvalidLevel()) {
          //G4cout << "Mismatch in angular momentum at E_level = " << E_level << "." << G4endl;
          //G4cout << "_angularMomentum = " << _angularMomentum << G4endl;
          //G4cout << "spin = " << spin << G4endl;
          //G4cout << "Setting invalid level at E_level =" << E_level << G4endl;
          //G4cout << G4endl;
          this->SetInvalidLevel();
        } else {
          if (_Verbose) {
            G4cout << "Spins agree: " << spin << G4endl;
          }
        }
        _Tau           = E_width;
        _halfLife      = T_half;
        _parity        = parity;
        _Ewidth_gamma  *= eV;
        _Ewidth_gamma0 *= eV;
        _Ewidth_p      *= eV;
        _Ewidth_n      *= eV;
        _Ewidth_alpha  *= eV;
      } // close (if E difference <= tol) check
    }
  } // close while loop

  if (!found_level) {
    if (false) {
      G4cout << "RefreshWidth: Could not find level " << _energy/keV << " keV. " << G4endl;
      G4cout << "  no such level in " << levels_filename << G4endl;
      G4cout << "Setting invalid level" << G4endl;
      G4cout << G4endl;
    }
    //G4cout << "Setting invalid level at E/MeV = " << _energy/MeV << G4endl;
    this->SetInvalidLevel();
  }
}


void G4NRFNuclearLevel::MakeWidth0() {
  // NRF-specific routine
  //
  // Calculates NRF gs transition partial width (_Tau0) of current level.
  // This routine enforces a policy for calculating _Tau0 based on
  // available ENSDF info.  Calls Calc_Tau0() to determine _Tau0
  // from value of electromagnetic decay width of level, _Ewidth_gamma.


  // Policy for determining ground-state partial width of level (_Tau0):
  //
  // If natural linewidth of level (_Tau) is unknown, don't permit
  // NRF excitation of state (set _Tau0 = 0).
  //
  // If natural linewidth is known and the gamma-decay width is explicitly
  // tabulated in ENSDF, then calculate gs width using gammay decay
  // probabilities (call Calc_Tau0()).
  //
  // If natural linewidth is known but the gamma-decay width is NOT
  // stated in ENSDF, then determine from the level width if it
  // appears to decay by gamma emission.  If so, then ASSUME that
  // the gamma-decay width equals the natural linewidth and
  // calculate the gs-width as above.  If the width of the state
  // does not appear to be compatible with EM transition, then
  // do not permit NRF excitation (set _Tau0 = 0).

  const G4double TAU_CUT_LO = 1.e-10 * eV;
  const G4double TAU_CUT_HI = 1.0 * eV;

  if (_Tau < 0.0) {   // natural linewidth unknown
    _Tau0 = 0.0;      // do not permit NRF excitation
  } else {    // _Tau is known
    if (_Ewidth_gamma > 0.0) {      // gamma decay width in ENSDF
      Calc_Tau0();                  // calculate gs width
    } else {                        // gamma decay width NOT in ENSDF
      if ((_Tau >= TAU_CUT_LO) &&   // check if level "looks" like
          (_Tau <= TAU_CUT_HI)) {   // it decays primarily by EM process
        _Ewidth_gamma = _Tau;       // If yes, assume gamma decay width = natural width
        Calc_Tau0();                // And calculate gs width
      } else {     // Level width incompatible with EM transition & nothing known
                   // about gamma width -- don't permit NRF excitation
        _Tau0 = 0.0;
        //  G4cout << "Level width of " << _Tau/eV << " eV falls outside NRF tau cut."
        //         << G4endl;
      }
    }
  }
}


void G4NRFNuclearLevel::Calc_Tau0() {
  // NRF-specific routine
  //
  // Precondition: Member variables _nGammas, _weights[], _Ewidth_gamma, _prob[],
  // and _totalCC[] must be previously initialized.
  //
  // Calculates gs gamma decay width (_Tau0) from total electromagnetic width using
  // relative gamma emission probabilities previously calculated and stored in the member
  // variable array _prob[].  N.B. The transition probabilities calculated in MakeProbabilities()
  // refer to total (gamma + conversion electron) emissions.  The _Tau0 calculated here
  // refers to gamma emission only.
  //
  // N.B. Function assumes that the variable _Ewidth_gamma (electromagnetic decay width of level)
  // has been previously set in MakeWidth0().
  //
  // Gamma transition to gs is identified via call to Identify_GS_Transition().
  //
  // If no gammas to gs occur from the level, set _Tau0 = 0.

  _Tau0 = 0.0;  // Default to zero (in case no gammas to gs)


  if (_nGammas > 0) {
    G4bool gamma_to_gs = Identify_GS_Transition();

    // N.B. Gamma "weights" can be zero, whereas
    // gamma "probs" are unit-normalized.

    if (gamma_to_gs) {
      if (_weights[_nGammas-1] > 0.0) {
        _Tau0 = _Ewidth_gamma * _prob[_nGammas-1]/(1+ _totalCC[_nGammas-1]);
      } else {
        _Tau0 = 0.0;
      }

      if (_Verbose) {
        G4cout << "G4NRFNuclearLevel::Calc_Tau0(): Calculated _Tau0 (eV) = "
               << _Tau0/eV << G4endl;
        G4cout << " from _Ewidth_gamma (eV): " << _Ewidth_gamma/eV
               << " Gamma probability: " << _prob[_nGammas-1]
               << " weight: " << _weights[_nGammas-1]
               << " Total CC: " << _totalCC[_nGammas-1]
               << " _nGammas: " << _nGammas << G4endl;
      }
    }
  }
}


G4bool G4NRFNuclearLevel::Identify_GS_Transition() {
  // NRF-specific routine
  //
  // Precondition: _nGammas > 0
  // Member variables _energy, _E_1stExcitedState, _nucleusA, _energies[]
  // must be previously initialized.
  //
  // Determines whether the level's highest-energy gamma (i.e.,
  //             E_gamma = _energies[_nGammas-1] )
  // corresponds to transition to the ground state.  The function calculates
  // the recoil-corrected gamma emission energy for (a) a transition to the gs,
  // and (b) a transition to the 1st excited state.  If the absolute value
  // of the energy difference (actual gamma energy - gamma emission energy)
  // is closer to (a), then the function returns true (means "transition
  // to gs").  Otherwise returns false (means "the level does not have a
  // gamma transition to gs").


  G4bool gamma_to_gs = false;  // Default: NOT a transition to gs

  // Use CLHEP::amu rather than the following:
  // const G4double amu = 931.5016 * MeV;
  // G4cout << "Using CLHEP::amu_c2 (MeV) = " << amu_c2/MeV << G4endl;

  G4double E_level = _energy;
  G4double E_1st   = _E_1stExcitedState;
  G4double M0      = _nucleusA * amu_c2;

  if (_Verbose) {
    G4cout << "!!!!! In Identify_GS_Transition().  E_level (MeV): "
     << E_level/MeV << G4endl;
    G4cout << "      E_1st (MeV): " << E_1st/MeV << G4endl;
  }

  // Note: Gammas are listed in order of increasing energy,
  // so gamma to gs (if it occurs) should always be the last in the
  // list.
  G4double E_gamma = _energies[_nGammas-1];


  // Nuclear recoil energy (non-relativistic approx.); see G4NRF class for factor of 2 explanation
  G4double E_recoil = E_level*E_level/(2.0 * M0);
  G4double E_gamma_to_gs = E_level - E_recoil;

  if (_Verbose) {
    G4cout << "   E_gamma_to_gs (keV): " << E_gamma_to_gs/keV << G4endl;
    G4cout << "   E_gamma (keV):       " << E_gamma/keV << G4endl;
  }

  G4double Ediff_gs = fabs(E_gamma - E_gamma_to_gs);

  if (_Verbose) {
    G4cout << "   Ediff_gs (keV):      " << Ediff_gs/keV << G4endl;
  }

  // Repeat calculation for 1st excited state
  G4double deltaE = E_level - E_1st;  // "zero'th order" E_gamma to 1st exc. state
  G4double M1 = M0 + E_1st;           // mass of 1st exc. state
  E_recoil = deltaE*deltaE/(2.0 * M1);
  G4double E_gamma_to_1st = deltaE - E_recoil;
  G4double Ediff_1st = fabs(E_gamma - E_gamma_to_1st);

  if (_Verbose) {
    G4cout << "   E_gamma_to_1st (keV): " << E_gamma_to_1st/keV << G4endl;
    G4cout << "   Ediff_1st (keV):      " << Ediff_1st/keV << G4endl;
  }

  // Is the actual gamma energy closer to gamma emission energy to gs or 1st exc.?
  if (Ediff_gs < Ediff_1st)
    gamma_to_gs = true;

  if (_Verbose) {
    G4cout << "   Returning Bool gamma_to_gs: " << gamma_to_gs << G4endl;
  }

  return gamma_to_gs;
}

const interpolating_function_p<G4double>* G4NRFNuclearLevel::GetCrossSectionTable() const {return _cross_sec_interp_func;}

void G4NRFNuclearLevel::SetCrossSectionTable(interpolating_function_p<G4double> *f) {_cross_sec_interp_func = f;}
