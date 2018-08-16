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
// Description:
//
// G4NRFNuclearLevelStore is a singleton used to contain all the G4NRFNuclearLevelManager objects.
//
// Updates by Jayson Vavrek (jvavrek@mit.edu):
//
// 0) Replaced the previous exceptionally slow look-up table system. Rather than looking up databases
//    at each cross section evaluation using several slow string comparisons and concatenations, use
//    a much faster integer key system that is created at initialization.
//
// -------------------------------------------------------------------

#include "G4NRFNuclearLevelStore.hh"
#include <sstream>

std::map<G4String, G4NRFNuclearLevelManager*> G4NRFNuclearLevelStore::theManagers;

// much faster lookup table system for getting the Managers based on an integer rather than a string comp
std::vector<G4String> G4NRFNuclearLevelStore::theKeys_fast(300*100, "");
std::vector<G4NRFNuclearLevelManager*> G4NRFNuclearLevelStore::theManagers_fast(300*100, NULL);

G4String G4NRFNuclearLevelStore::dirName("");

G4NRFNuclearLevelStore* G4NRFNuclearLevelStore::GetInstance() {
  static G4NRFNuclearLevelStore theInstance;
  return &theInstance;
}

G4NRFNuclearLevelStore::G4NRFNuclearLevelStore() {
  char* env = getenv("G4NRFGAMMADATA");
  if (!env) {
    G4cout << "G4NRFNuclearLevelStore: please set the G4NRFGAMMADATA environment variable"
           << G4endl;
    G4cout << "Aborting." << G4endl;
    exit(51);
  } else {
    dirName = env;
    dirName += '/';
  }

  for (int a = 1; a <= 300; a++) {
    for (int z = 1; z <= 100; z++) {
      theKeys_fast[GetKeyIndex(z, a)] = this->GenerateKey(z, a);
    }
  }
}


G4NRFNuclearLevelStore::~G4NRFNuclearLevelStore() {
  std::map<G4String, G4NRFNuclearLevelManager*>::iterator i;
  for (i = theManagers.begin(); i != theManagers.end(); ++i) {
    if ( (*i).second ) delete (*i).second;
  }
}

// this was the time-limiting function call in the original code -- too many slow string ops per step
G4String G4NRFNuclearLevelStore::GenerateKey(const G4int Z, const G4int A) {
  std::ostringstream streamName;
  streamName << 'z' << Z << ".a" << A;
  G4String name(streamName.str());
  return name;
}


G4int G4NRFNuclearLevelStore::GetKeyIndex(const G4int Z, const G4int A) {
  return 100*(A-1) + (Z-1);
}


G4NRFNuclearLevelManager* G4NRFNuclearLevelStore::GetManager(const G4int Z, const G4int A, G4bool standalone) {
  G4NRFNuclearLevelManager * result = 0;
  if (A < 1 || Z < 1 || A < Z) {
    G4cerr << "G4NRFNuclearLevelStore::GetManager: Wrong values Z = " << Z << " A = " << A << '\n';
    return result;
  }
  // Generate the key = filename
  G4int index = GetKeyIndex(Z, A);
  G4String key(theKeys_fast[index]);

  // Check if the manager exists. If it does, return it. Otherwise, create it.
  G4bool managerExists = (theManagers_fast[index] != NULL ? true : false);
  if (!managerExists) {
    result = new G4NRFNuclearLevelManager();
    result->SetNucleus(Z, A, dirName + key, standalone);
    theManagers_fast[index] = result;
  } else {
    result = theManagers_fast[index];
  }

  return result;
}
