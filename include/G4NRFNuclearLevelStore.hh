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

#ifndef G4NRFNuclearLevelStore_hh
#define G4NRFNuclearLevelStore_hh 1

#include <map>
#include <vector>

#include "G4NRFNuclearLevelManager.hh"

class G4NRFNuclearLevel;

class G4NRFNuclearLevelStore {
 private:
  G4NRFNuclearLevelStore();

 public:
  static G4NRFNuclearLevelStore* GetInstance();

  G4NRFNuclearLevelManager * GetManager(const G4int Z, const G4int A, G4bool standalone = false);

  ~G4NRFNuclearLevelStore();

 private:
  G4String GenerateKey(const G4int Z, const G4int A);

  G4int GetKeyIndex(const G4int Z, const G4int A);
  static std::vector<G4String> theKeys_fast;
  static std::vector<G4NRFNuclearLevelManager*> theManagers_fast;

  static std::map<G4String, G4NRFNuclearLevelManager*> theManagers;
  static G4String dirName;
};

#endif
