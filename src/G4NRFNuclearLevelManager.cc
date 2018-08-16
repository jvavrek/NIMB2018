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
//      File name:     G4NRFNuclearLevelManager
//
//      Author:        David Jordan (david.jordan@pnl.gov)
//      Creation date: October 2006
//
//      Adapted from   G4NuclearLevelManager
//      CERN, Geneva, Switzerland
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//      Creation date: 24 October 1998
//
//      Modification history (G4NuclearLevelManager):
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//        02 May 2003,   Vladimir Ivanchenko remove public copy constructor
//        14 Jan 2015,   Jayson Vavrek (jvavrek@mit.edu)
//              Moderate to major bug fixes plus enhanced documentation.
//
//
// -------------------------------------------------------------------

// -------------------------------------------------------------------
//
// Description:
//
// G4NRFNuclearLevelManager is the class used to import nuclear level data
// from databases and build instantiations of G4NRFNuclearLevel.
//
//
// Updates by Jayson Vavrek:
//
// 0) Modernized (2015) C++ standards: #include's, using's; also
//    G4SystemOfUnits.hh and G4PhysicalConstants.hh
//
// 1) Levels with empirically-unknown spins or badly-mismatched energies are
//    erased via delete_bad_levels(), eliminating entirely the angular momentum
//    and mismatch errors present earlier. Note that this may delete useful levels
//    and the user must edit databases by hand in response.
//
// 2) Added standalone NRF gamma database generator functionality. The standalone
//    bool gets passed from G4NRF.cc, and if true, disables some error checking
//    (which may be redundant with new level-matching scheme anyway) so that the
//    database can be compiled without exit() errors. The PrintAllTabular() method
//    uses G4NRFNuclearLevel::PrintAllTabular() to print everything to an ofstream
//    file that also gets set in G4NRF.cc.
//
// -------------------------------------------------------------------

#include "G4NRFNuclearLevelManager.hh"

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <typeinfo>

#include "G4NRFNuclearLevelStore.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4HadTmpUtil.hh"
#include "G4HadronicException.hh"

G4NRFNuclearLevelManager::G4NRFNuclearLevelManager(G4bool Verbose):
  _nucleusA(0), _nucleusZ(0), _fileName(""), _validity(false),
  _levels(0), _levelEnergy(0), _gammaEnergy(0), _probability(0),
  _Verbose(Verbose) { }

G4NRFNuclearLevelManager::G4NRFNuclearLevelManager(const G4int Z, const G4int A, const G4String& filename, G4bool Verbose) :
  _nucleusA(A), _nucleusZ(Z), _fileName(filename), _Verbose(Verbose) {
  if (A <= 0 || Z <= 0 || Z > A)
    throw G4HadronicException(__FILE__, __LINE__, "==== G4NRFNuclearLevelManager ==== (Z, A) < 0, or Z > A");

  _levels = 0;

  if (_Verbose) {
    G4cout << "G4NRFNuclearLevelManager: Initializing (Z, A) = ("
     << Z << ","
     << A << ")" << G4endl;
  }

  ReadGroundStateProperties(); // never gets called

  MakeLevels(); // never happens because this particular constructor never gets called
}

G4NRFNuclearLevelManager::~G4NRFNuclearLevelManager() {
  if (_levels) {
    if (_levels->size() > 0) {
      std::for_each(_levels->begin(), _levels->end(), DeleteLevel());
      _levels->clear();
    }

    delete _levels;
    _levels = 0;
  }
}

void G4NRFNuclearLevelManager::SetNucleus(const G4int Z, const G4int A, const G4String& filename, G4bool standalone) {
  if (_nucleusZ != Z || _nucleusA != A) {
    _nucleusA = A;
    _nucleusZ = Z;
    _fileName = filename;

    if (_Verbose) {
      G4cout << "G4NRFNuclearLevelManager::SetNucleus Initializing (Z, A) = ("
             << Z << ","
             << A << ")" << G4endl;

      G4cout << "G4NRFNuclearLevelManager::SetNucleus Calling ReadGroundStateProperties"
             << G4endl;
    }

    if (_Verbose) {
      G4cout << "G4NRFNuclearLevelManager::SetNucleus Calling MakeLevels"
             << G4endl;
    }

    MakeLevels();
    ReadGroundStateProperties(standalone); // this is the only place ReadGroundStateProperties ever gets called
    ReadTDebyeData(standalone);
    const G4double Teff = CalcTeff(_TDebye);
    SetTeff(Teff);

    if (_Verbose) {
      G4cout << "G4NRFNuclearLevelManager::SetNucleus Back from MakeLevels" << G4endl;
    }
  }
}

void G4NRFNuclearLevelManager::ReadGroundStateProperties(G4bool standalone) {
  // NRF-specific routine
  //
  // Reads spin and parity of ground state from file
  // $G4NRFGAMMADATA/ground_state_properties.dat

  char* env = getenv("G4NRFGAMMADATA");
  if (!env) {
    G4cout << "G4NRFNuclearLevelManager: please set the G4NRFGAMMADATA environment variable"
           << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(2);
  }

  G4String dir = G4String(env) + '/';

  G4String gs_filename = dir + G4String("ground_state_properties.dat");

  std::ifstream gsFile(gs_filename, std::ios::in);
  if (!gsFile) {
    G4cout << "Error in G4NRFNuclearLevelManager::ReadGroundStateProperties." << G4endl;
    G4cout << "Could not open file ground_state_properties.dat" << G4endl;
    G4cout << "Expected to find this file in directory" << G4endl;
    G4cout << dir << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(3);
  }

  G4bool found_gs = false;

  char line[80];

  G4int Z, A;
  G4double spin, parity;

  while ((!found_gs) && (!gsFile.eof())) {
    gsFile.getline(line, 80);
    std::istringstream line_str(line);

    line_str >> Z >> A >> spin >> parity;

    if ((Z == _nucleusZ) && (A == _nucleusA)) {
      found_gs = true;
      if (spin >= 0.0) {
        _gsAngularMomentum = spin;

        if (_Verbose) {
          G4cout << "G4NRFNuclearLevelManager::ReadGroundStateProperties -- Found valid g.s. spin = "
                 << spin << G4endl;
        }
      } else {
        G4cout << "WARNING from G4NRFNuclearLevelManager::ReadGroundStateProperties."
               << G4endl;
        G4cout << "Ground state spin is unknown, setting to 0..."
               << G4endl;
        _gsAngularMomentum = 0.0;
      }

      if (parity != 0.0) {
        _gsParity = parity;
        if (_Verbose) {
          G4cout << "G4NRFNuclearLevelManager::ReadGroundStateProperties -- Found valid g.s. parity = "
                 << parity << G4endl;
        }
      } else {
        G4cout << "WARNING from G4NRFNuclearLevelManager::ReadGroundStateProperties."
               << G4endl;
        G4cout << "Ground state parity is unknown, setting to +"
               << G4endl;
        _gsParity = 1.0;
      }
    }
  }

  gsFile.close();

  if (!found_gs) {
    G4cout << "ERROR in G4NRFNuclearLevelManager::ReadGroundStateProperties." << G4endl;
    G4cout << "Could not find g.s. for Z = " << _nucleusZ << " A = " << _nucleusA << G4endl;
    if (!standalone) G4cout << "Aborting." << G4endl;
    if (!standalone) exit(4);
    if (standalone) G4cout << "Skipping ground state for standalone code" << G4endl;
  }
}

void G4NRFNuclearLevelManager::ReadTDebyeData(G4bool standalone) {
  char* env = getenv("G4NRFGAMMADATA");
  if (!env) {
    G4cout << "G4NRFNuclearLevelManager: please set the G4NRFGAMMADATA environment variable" << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(5);
  }

  G4String dir = G4String(env) + '/';

  G4String TD_filename = dir + G4String("TDebye_data.dat");

  std::ifstream TDFile(TD_filename, std::ios::in);
  if (!TDFile) {
    G4cout << "Error in G4NRFNuclearLevelManager::ReadGroundStateProperties." << G4endl;
    G4cout << "Could not open file TDebye_data.dat" << G4endl;
    G4cout << "Expected to find this file in directory" << G4endl;
    G4cout << dir << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(6);
  }

  G4bool found_TD = false;
  char line[80];
  G4double TDA, TDK; // T_Debye values from Ashcroft/Mermin and Kittel, respectively
  G4int Z;           // Z of nucleus
  G4bool LA, LK;     // whether the T_Debye value is a low-temperature determination
  G4double TDebye_tmp = 0.0; // initialize temporary TDebye value

  while ((!found_TD) && (!TDFile.eof())) {
    TDFile.getline(line, 80);
    std::istringstream line_str(line);

    line_str >> Z >> TDA >> LA >> TDK >> LK;

    if (Z == _nucleusZ) {
      if        (TDA == 0 && TDK == 0) {
        TDebye_tmp = 0.0;
      } else if (TDA == 0 && TDK != 0) {
        TDebye_tmp = TDK;
      } else if (TDA != 0 && TDK == 0) {
        TDebye_tmp = TDA;
      } else {
        if      ( LA && !LK) {
          TDebye_tmp = TDK;
        } else if (!LA &&  LK) {
          TDebye_tmp = TDA;
        } else {
          TDebye_tmp = (TDA + TDK)/2.0;
        }
      }

      found_TD = true;
    }
  }

  TDFile.close();

  if (!found_TD) {
    G4cout << "ERROR in G4NRFNuclearLevelManager::ReadTDebyeData." << G4endl;
    G4cout << "Could not find T_Debye for Z = " << Z << G4endl;
    if (!standalone) G4cout << "Aborting." << G4endl;
    if (!standalone) exit(8);
    if (standalone) G4cout << "Skipping ground state for standalone code" << G4endl;
  } else {
    SetTDebye(TDebye_tmp*kelvin);
  }
}

void G4NRFNuclearLevelManager::SetTDebye(G4double TDebye) {
  _TDebye = TDebye;
}

G4double G4NRFNuclearLevelManager::GetTDebye() const {
  return _TDebye;
}

G4double G4NRFNuclearLevelManager::GetGroundStateSpin() const {
  return _gsAngularMomentum;
}

G4double G4NRFNuclearLevelManager::GetGroundStateParity() const {
  return _gsParity;
}

G4bool G4NRFNuclearLevelManager::IsValid() const {
  return _validity;
}

G4int G4NRFNuclearLevelManager::NumberOfLevels() const {
  G4int n = 0;
  if (_levels != 0) n = _levels->size();
  return n;
}

const G4NRFPtrLevelVector* G4NRFNuclearLevelManager::GetLevels() const {
  return _levels;
}


G4NRFNuclearLevel* G4NRFNuclearLevelManager::
NearestLevelRecoilAbsorb(const G4double GammaEnergy, const G4double eDiffMax) {
  // NRF-specific routine
  //
  // Find closest nuclear level, correcting for nuclear recoil
  // upon absorption of gamma.  Gamma energy must equal sum of
  // internal excitation energy (i.e. level energy) plus
  // kinetic energy imparted to nucleus.  Nucleus is assumed
  // to be in g.s. initially.
  //
  // Nuclear mass used for (non-relativistic) recoil calculation
  // neglects mass defect.

  const G4double M = _nucleusA * amu_c2;

  const G4double E_recoil = GammaEnergy*GammaEnergy/(2.0 * M); // see G4NRF class for factor of 2 explanation

  G4double E_internal_excite = GammaEnergy - E_recoil;

  if (_Verbose) {
    G4cout << "G4NRFNuclearLevelManager::NearestLevelRecoilAbsorb -- " << G4endl;
    G4cout << "GammaEnergy (MeV): "       << GammaEnergy/MeV       << G4endl;
    G4cout << "Atomic weight:     "       << _nucleusA             << G4endl;
    G4cout << "E_recoil (keV):    "       << E_recoil/keV          << G4endl;
    G4cout << "E_internal_excite (MeV): " << E_internal_excite/MeV << G4endl;
  }

  return NearestLevel(E_internal_excite, eDiffMax);
}

G4NRFNuclearLevel* G4NRFNuclearLevelManager::
NearestLevelRecoilEmit(const G4double LevelEnergy, const G4double GammaEnergy, const G4double eDiffMax) {
  // NRF-specific routine
  //
  // Find closest nuclear level, correcting for nuclear recoil
  // upon emission of gamma from level of energy LevelEnergy.
  // Change in internal excitation must equal sum of emitted
  // gamma energy plus recoil energy imparted to nucleus.
  //
  // Nuclear mass used for (non-relativistic) recoil calculation
  // neglects mass defect.

  const G4double M = _nucleusA * amu_c2;

  const G4double E_recoil = GammaEnergy*GammaEnergy/(2.0 * M); // see G4NRF class for factor of 2 explanation

  G4double E_internal_DeExcite = GammaEnergy + E_recoil;

  G4double Target_LevelEnergy = LevelEnergy - E_internal_DeExcite;

  return NearestLevel(Target_LevelEnergy, eDiffMax);
}


G4NRFNuclearLevel* G4NRFNuclearLevelManager::
NearestLevel(const G4double energy, const G4double eDiffMax) {
  G4int iNear = -1;
  G4double diff = 9999. * GeV;

  if (_levels != 0) {
    G4double eLast = 0.0;
    G4double eDiffLast = diff;
    unsigned int i = 0;
    const int numLevels = _levels->size(); // this saves almost a full 1% CPU
    for (i = 0; i < numLevels; ++i) {
      G4double e = _levels->operator[](i)->Energy();

      if (e <= eLast) {G4cout << "error: lower energy!" << G4endl; exit(18);} // if the Energy vals aren't ordered, throw error

      G4double eDiff = std::abs(e - energy);
      if (eDiff < diff && eDiff <= eDiffMax) {
        diff = eDiff;
        iNear = i;
      }

      // the Energy() vals appear to be ordered, so if eDiff starts increasing, terminate loop
      if (eDiff > eDiffLast && i != 0) break;
      eDiffLast = eDiff;
    }
  }

  if (_levels != 0 && iNear >= 0 && iNear < static_cast<G4int>(_levels->size()) ) {
    return _levels->operator[](iNear);
  } else {
    return 0;
  }
}


G4double G4NRFNuclearLevelManager::MinLevelEnergy() const {
  G4double eMin = 9999.*GeV;
  if (_levels != 0) {
    if (_levels->size() > 0) eMin = _levels->front()->Energy();
  }
  return eMin;
}


G4double G4NRFNuclearLevelManager::MaxLevelEnergy() const {
  G4double eMax = 0.;
  if (_levels != 0) {
    if (_levels->size() > 0) eMax = _levels->back()->Energy();
  }
  return eMax;
}


const G4NRFNuclearLevel* G4NRFNuclearLevelManager::HighestLevel() const {
  if (_levels!= 0 && _levels->size() > 0)
    return _levels->front();
  else
    return 0;
}


const G4NRFNuclearLevel* G4NRFNuclearLevelManager::LowestLevel() const {
  if (_levels != 0 && _levels->size() > 0)
    return _levels->back();
  else
    return 0;
}


G4bool G4NRFNuclearLevelManager::Read(std::ifstream& dataFile) {
  const G4double minProbability = 0.001;

  G4bool result = true;

  if (dataFile >> _levelEnergy) {
    dataFile >> _gammaEnergy >> _probability >> _polarity >> _halfLife
             >> _angularMomentum  >> _totalCC >> _kCC >> _l1CC >> _l2CC
             >> _l3CC >> _m1CC >> _m2CC >> _m3CC >> _m4CC >> _m5CC
             >> _nPlusCC;

    _levelEnergy *= keV;
    _gammaEnergy *= keV;
    _halfLife *= second;

    // Comment appearing in original G4NuclearLevelManager.cc:
    // The following adjustment is needed to take care of anomalies in
    // data files, where some transitions show up with relative probability
    // zero

    // Note added (DJ 28-Aug-06): The "anomalies" in the probabilities
    // evidently correspond to forbidden gamma transitions.  Disable
    // the provision for setting threshold "minProbability" in order to
    // eliminate NRF transitions to these forbidden levels.

    //  if (_probability < minProbability) _probability = minProbability;

    // the folowwing is to convert icc probability to accumulative ones
    _l1CC += _kCC;
    _l2CC += _l1CC;
    _l3CC += _l2CC;
    _m1CC += _l3CC;
    _m2CC += _m1CC;
    _m3CC += _m2CC;
    _m4CC += _m3CC;
    _m5CC += _m4CC;
    _nPlusCC += _m5CC;
    if (_nPlusCC != 0) {
      _kCC /= _nPlusCC;
      _l1CC /= _nPlusCC;
      _l2CC /= _nPlusCC;
      _l3CC /= _nPlusCC;
      _m1CC /= _nPlusCC;
      _m2CC /= _nPlusCC;
      _m3CC /= _nPlusCC;
      _m4CC /= _nPlusCC;
      _m5CC /= _nPlusCC;
      _nPlusCC /= _nPlusCC;
    } else {
      _kCC = 1;
      _l1CC = 1;
      _l2CC = 1;
      _l3CC = 1;
      _m1CC = 1;
      _m2CC = 1;
      _m3CC = 1;
      _m4CC = 1;
      _m5CC = 1;
      _nPlusCC = 1;
    }
  } else {
    result = false;
  }

  return result;
}


void G4NRFNuclearLevelManager::MakeLevels() {
  _validity = false;

  std::ifstream inFile(_fileName, std::ios::in);
  if (!inFile) {
    if (_nucleusZ > 10) G4cout << " G4NRFNuclearLevelManager: nuclide ("
                               << _nucleusZ << "," << _nucleusA
                               << ") does not have a gamma levels file" << G4endl;
    return;
  }

  if (_levels != 0) {
    if (_levels->size() > 0) {
      std::vector<G4NRFNuclearLevel*>::iterator pos;
      for (pos = _levels->begin(); pos != _levels->end(); pos++)
        if (*pos) delete *pos;
      _levels->clear();
    }
    delete _levels;
  } else {
    _validity = true;
  }

  _levels = new G4NRFPtrLevelVector;

  std::vector<G4double> eLevel;
  std::vector<G4double> eGamma;
  std::vector<G4double> wGamma;
  std::vector<G4double> pGamma; // polarity
  std::vector<G4double> hLevel; // half life
  std::vector<G4double> aLevel; // angular momentum
  std::vector<G4double> kConve; // internal convertion coefficiencies
  std::vector<G4double> l1Conve;
  std::vector<G4double> l2Conve;
  std::vector<G4double> l3Conve;
  std::vector<G4double> m1Conve;
  std::vector<G4double> m2Conve;
  std::vector<G4double> m3Conve;
  std::vector<G4double> m4Conve;
  std::vector<G4double> m5Conve;
  std::vector<G4double> npConve;
  std::vector<G4double> toConve;


  while (Read(inFile)) {
    eLevel.push_back(_levelEnergy);
    eGamma.push_back(_gammaEnergy);
    wGamma.push_back(_probability);
    pGamma.push_back(_polarity);
    hLevel.push_back(_halfLife);
    aLevel.push_back(_angularMomentum);
    kConve.push_back(_kCC);
    l1Conve.push_back(_l1CC);
    l2Conve.push_back(_l2CC);
    l3Conve.push_back(_l3CC);
    m1Conve.push_back(_m1CC);
    m2Conve.push_back(_m2CC);
    m3Conve.push_back(_m3CC);
    m4Conve.push_back(_m4CC);
    m5Conve.push_back(_m5CC);
    npConve.push_back(_nPlusCC);
    toConve.push_back(_totalCC);
  }

  // ---- MGP ---- Don't forget to close the file
  inFile.close();

  G4int nData = eLevel.size();

  if (_Verbose) {
    G4cout << " ==== MakeLevels ===== " << nData << " data read " << G4endl;
  }

  G4double thisLevelEnergy = eLevel[0];
  G4double E_1stExcitedState = eLevel[0];  // Added DJ 8/29/06
  G4double thisLevelHalfLife = 0.;
  G4double thisLevelAngMom = 0.;
  std::vector<G4double> thisLevelEnergies;
  std::vector<G4double> thisLevelWeights;
  std::vector<G4double> thisLevelPolarities;
  std::vector<G4double> thisLevelkCC;
  std::vector<G4double> thisLevell1CC;
  std::vector<G4double> thisLevell2CC;
  std::vector<G4double> thisLevell3CC;
  std::vector<G4double> thisLevelm1CC;
  std::vector<G4double> thisLevelm2CC;
  std::vector<G4double> thisLevelm3CC;
  std::vector<G4double> thisLevelm4CC;
  std::vector<G4double> thisLevelm5CC;
  std::vector<G4double> thisLevelnpCC;
  std::vector<G4double> thisLeveltoCC;

  G4double e = -1.;
  G4int i;
  G4int jlevel = 0; // start levels at 1=1st excited state

  for (i = 0; i < nData; i++) {
    e = eLevel[i];
    if (e != thisLevelEnergy) {
      //    G4cout << "Making a new level... " << e << G4endl;
      //     << thisLevelEnergies.entries() << " "
      //     << thisLevelWeights.entries() << G4endl;

      // DJ 8/29/06  Added E_1stExcitedState to G4NRFNuclearLevel
      // constructor to aid in identifying gamma transitions
      // to ground state (distinction between gamma to gs and
      // 1st exc. state may not be clear otherwise for heavy
      // nuclei)

      jlevel++;

      // this is where most of the heavy lifting should get done?
      // ideally, skip this step altogether if the spin = -999
      G4NRFNuclearLevel* newLevel = new G4NRFNuclearLevel(jlevel,
              _nucleusZ,
              _nucleusA,
              E_1stExcitedState,
              thisLevelEnergy,
              thisLevelHalfLife,
              thisLevelAngMom,
              thisLevelEnergies,
              thisLevelWeights,
              thisLevelPolarities,
              thisLevelkCC,
              thisLevell1CC,
              thisLevell2CC,
              thisLevell3CC,
              thisLevelm1CC,
              thisLevelm2CC,
              thisLevelm3CC,
              thisLevelm4CC,
              thisLevelm5CC,
              thisLevelnpCC,
              thisLeveltoCC);

      _levels->push_back(newLevel);
      // Reset data vectors
      thisLevelEnergies.clear();
      thisLevelWeights.clear();
      thisLevelPolarities.clear();
      thisLevelkCC.clear();
      thisLevell1CC.clear();
      thisLevell2CC.clear();
      thisLevell3CC.clear();
      thisLevelm1CC.clear();
      thisLevelm2CC.clear();
      thisLevelm3CC.clear();
      thisLevelm4CC.clear();
      thisLevelm5CC.clear();
      thisLevelnpCC.clear();
      thisLeveltoCC.clear();
      thisLevelEnergy = e;
    }

    // Append current data
    thisLevelEnergies.push_back(eGamma[i]);
    thisLevelWeights.push_back(wGamma[i]);
    thisLevelPolarities.push_back(pGamma[i]);
    thisLevelkCC.push_back(kConve[i]);
    thisLevell1CC.push_back(l1Conve[i]);
    thisLevell2CC.push_back(l2Conve[i]);
    thisLevell3CC.push_back(l3Conve[i]);
    thisLevelm1CC.push_back(m1Conve[i]);
    thisLevelm2CC.push_back(m2Conve[i]);
    thisLevelm3CC.push_back(m3Conve[i]);
    thisLevelm4CC.push_back(m4Conve[i]);
    thisLevelm5CC.push_back(m5Conve[i]);
    thisLevelnpCC.push_back(npConve[i]);
    thisLeveltoCC.push_back(toConve[i]);
    thisLevelHalfLife = hLevel[i];
    thisLevelAngMom = aLevel[i];
  }

  // Make last level
  if (e > 0.) {
    jlevel++;

    G4NRFNuclearLevel* newLevel = new G4NRFNuclearLevel(jlevel,
                  _nucleusZ,
                  _nucleusA,
                  E_1stExcitedState,
                  e,
                  thisLevelHalfLife,
                  thisLevelAngMom,
                  thisLevelEnergies,
                  thisLevelWeights,
                  thisLevelPolarities,
                  thisLevelkCC,
                  thisLevell1CC,
                  thisLevell2CC,
                  thisLevell3CC,
                  thisLevelm1CC,
                  thisLevelm2CC,
                  thisLevelm3CC,
                  thisLevelm4CC,
                  thisLevelm5CC,
                  thisLevelnpCC,
                  thisLeveltoCC);

    _levels->push_back(newLevel);
  }

  // check here if the level has -999 spin or angularMomentum; if so, erase() it
  delete_bad_levels();

  G4PtrSort<G4NRFNuclearLevel>(_levels);

  //PrintLevelEnergies();

  return;
}


void G4NRFNuclearLevelManager::PrintAll() {
  G4int nLevels = 0;
  if (_levels != 0) nLevels = _levels->size();

  G4cout << " ==== G4NRFNuclearLevelManager ==== (" << _nucleusZ << ", " << _nucleusA
         << ") has " << nLevels << " levels" << G4endl
         << "Highest level is at energy " << MaxLevelEnergy() << " MeV "
         << G4endl << "Lowest level is at energy " << MinLevelEnergy()
         << " MeV " << G4endl;
  G4cout << "Ground state spin & parity: " << _gsAngularMomentum
         << "  " << _gsParity << G4endl;

  G4int i = 0;
  for (i = 0; i < nLevels; i++) {
    _levels->operator[](i)->PrintAll();
  }
}


// routine used for printing NRF gamma info to standalone database
void G4NRFNuclearLevelManager::PrintAllTabular(ofstream& file) {
  G4int nLevels = 0;
  if (_levels != 0) nLevels = _levels->size();

  for (int i = 0; i < nLevels; i++) {
    double TDebye = G4NRFNuclearLevelStore::GetInstance()->GetManager(_nucleusZ, _nucleusA)->GetTDebye();
    _levels->operator[](i)->PrintAllTabular(file, _gsAngularMomentum, TDebye);
  }
}

void G4NRFNuclearLevelManager::PrintLevelEnergies() {
  G4cout << "Level energies [MeV] for Z = " << _nucleusZ << ", A = " << _nucleusA << G4endl;

  const int n = _levels->size();
  for (int i = 0; i < n; i++) {
    G4cout << _levels->operator[](i)->Energy()/MeV << G4endl;
  }
}


G4NRFNuclearLevelManager::G4NRFNuclearLevelManager(const G4NRFNuclearLevelManager &right) {
  _levelEnergy = right._levelEnergy;
  _gammaEnergy = right._gammaEnergy;
  _probability = right._probability;
  _polarity = right._polarity;
  _halfLife = right._halfLife;
  _angularMomentum = right._angularMomentum;
  _kCC = right._kCC;
  _l1CC = right._l1CC;
  _l2CC = right._l2CC;
  _l3CC = right._l3CC;
  _m1CC = right._m1CC;
  _m2CC = right._m2CC;
  _m3CC = right._m3CC;
  _m4CC = right._m4CC;
  _m5CC = right._m5CC;
  _nPlusCC = right._nPlusCC;
  _totalCC = right._totalCC;
  _nucleusA = right._nucleusA;
  _nucleusZ = right._nucleusZ;
  _fileName = right._fileName;
  _validity = right._validity;

  if (right._levels != 0) {
    _levels = new G4NRFPtrLevelVector;
    G4int n = right._levels->size();
    G4int i;
    for (i = 0; i < n; i++)
      _levels->push_back(new G4NRFNuclearLevel(*(right._levels->operator[](i))));

    G4PtrSort<G4NRFNuclearLevel>(_levels);
  } else {
    _levels = 0;
  }
}

void G4NRFNuclearLevelManager::delete_bad_levels() {
  //G4cout << "Calling delete_bad_levels() on Z = " << _nucleusZ << ", A = " << _nucleusA << G4endl;
  //G4cout << "Initial # of levels = " << _levels->size() << G4endl;
  G4int badcounter = 0;
  for (int i = 0; i < _levels->size(); i++) {
    G4NRFNuclearLevel *checkLevel = _levels->operator[](i);
    G4bool invalidLevel = checkLevel->GetInvalidLevel(); // check whether the level is invalid
    if (invalidLevel) {
      //G4cout << "Deleting bad level " << i << " at " << checkLevel->Energy()/MeV << " MeV." << G4endl;
      _levels->erase(_levels->begin()+i); // if so, delete it
      //G4cout << "New levels size = " << _levels->size() << G4endl;
      ++badcounter;
    }
  }
  //G4cout << "# of bad levels = " << badcounter << G4endl;
  //G4cout << "New # of levels = " << _levels->size() << G4endl;
  //G4cout << G4endl;

  //G4cout << "New list of levels = " << G4endl;
  //for (int i = 0; i < _levels->size(); i++)
    //G4cout << "Level energy = " << (_levels->operator[](i))->Energy()/MeV << " MeV." << G4endl;
}


G4double G4NRFNuclearLevelManager::TeffIntegrand(G4double t) const {
  return t*t*t * (1.0/(G4Exp(t) - 1.0) + 0.5);
}


G4double G4NRFNuclearLevelManager::CalcTeff(G4double TDebye, G4double Treal) {
  G4Integrator<const G4NRFNuclearLevelManager, G4double(G4NRFNuclearLevelManager::*)(G4double) const> integrator;

  G4double B = TDebye/Treal;

  const G4double theLowerLimit = 1.e-4;
  const G4double theUpperLimit = B;
  const G4int    nIterations   = 100;

  G4double sum = integrator.Simpson(this, &G4NRFNuclearLevelManager::TeffIntegrand,
                                    theLowerLimit, theUpperLimit, nIterations);

  G4double ratio = (3.0 / pow(B, 3)) * sum;

  G4double fTeff = Treal * ratio;

  return fTeff;
}


void G4NRFNuclearLevelManager::SetTeff(G4double Teff) { _Teff = Teff; }


G4double G4NRFNuclearLevelManager::GetTeff() { return _Teff; }




