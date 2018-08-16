/////////////////////////////////////////////////////////////////////////////////
//
// name: rootStorageManager.hh
// date: 19 Sep 14
// auth: Zach Hartwig
// mail: hartwig@psfc.mit.edu
//
// desc: The rootStorageManager class is designed to handle all
//       ROOT-related task: creation of ROOT objects (TFile's,
//       TTree's, etc), file name assignments, and writing data to
//       disk. The class is constructed as a meyer's singleton, and,
//       as such, is accessible from anywhere in the code where it is
//       important to store simulation data. The user should add
//       his/her own methods to the class to provide for data readout
//       into the appropriate ROOT structures (TTrees, TNtuples, etc).
//       Note the ROOT file naming, file initialization, and
//       file/object writing is handled via the associated
//       rootStorageManagerMessgenger class.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef rootStorageManager_hh
#define rootStorageManager_hh 1

#include <vector>
#include <string>
using std::string;

#include "G4Event.hh"
#include "G4Run.hh"

#include "TFile.h"
#include "TTree.h"

#include "rootStorageManagerMessenger.hh"

class rootStorageManager {
 public:
  rootStorageManager(int iind = 0);
  ~rootStorageManager();

  static rootStorageManager *GetInstance();

  void CreateROOTObjects();
  void WriteROOTObjects(G4bool emergencyWrite = false);
  void FillTree(G4double eD, G4double w, G4double eB, G4double eG);
  void FillTreeFoil(G4double e, G4double ei, G4double a, G4double w, G4bool n, G4int ci, G4bool rs);
  void FillTreeFoilInc(G4double e, G4double w);
  void FillTreePlateInc(G4double e, G4double w);
  void GenerateFileNames();

  void SetFileName(G4String FN) {ROOTFileName = FN;}
  G4String GetFileName() {return ROOTFileName;}

  G4bool CheckForROOTObjects() {return ROOTObjectsExist;}

 private:
  static rootStorageManager *theRootStorageManager;
  int find;

  G4String ROOTFileName;

  rootStorageManagerMessenger *theMessenger;

  TFile *ROOTFile;
  TTree *ROOTTree;
  TTree *ROOTTreeFoil;
  TTree *ROOTTreeFoilInc;
  TTree *ROOTTreePlateInc;
  G4bool ROOTObjectsExist;

  // variables written to ROOTTree
  G4double energyDeposition;
  G4double weight;
  G4double beamEnergy;
  G4double gammaEnergy;

  // and to ROOTTreeFoil
  G4double energy;
  G4double energyInitial;
  G4double angle;
  G4double weightF;
  G4bool isNRF;
  G4int creatorProcessIndex;
  G4bool resonanceSample;

  G4double foilEnergyInc, foilWeightInc;
  G4double plateEnergyInc, plateWeightInc;
};

#endif
