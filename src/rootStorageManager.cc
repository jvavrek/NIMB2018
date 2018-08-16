/////////////////////////////////////////////////////////////////////////////////
//
// name: rootStorageManager.cc
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


#include "rootStorageManager.hh"

// C++
#include <sstream>

// Geant4
#include "G4SDManager.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// ROOT
#include "TChain.h"

// ZK
#include "rootStorageManagerMessenger.hh"


rootStorageManager *rootStorageManager::theRootStorageManager = 0;


rootStorageManager *rootStorageManager::GetInstance()
{ return theRootStorageManager; }


rootStorageManager::rootStorageManager(int iind)
  : find(iind),
    ROOTFileName(""), ROOTFile(new TFile), ROOTTree(new TTree), ROOTTreeFoil(new TTree),
    ROOTTreeFoilInc(new TTree), ROOTTreePlateInc(new TTree),
    ROOTObjectsExist(false),
    energyDeposition(0.), weight(0.), beamEnergy(0.) {
  if (theRootStorageManager)
    G4Exception("rootStorageManager::rootStorageManager()",
      "ZKBremRootException001",
      FatalException,
      "The rootStorageManager was constructed twice!");
  else
    theRootStorageManager = this;

  theMessenger = new rootStorageManagerMessenger(this);
}


rootStorageManager::~rootStorageManager() {
  delete ROOTFile;
  delete theMessenger;
}


// Method that creates all the required ROOT objects before a G4 run
// commences. Note that this method should be called from the G4 cmd
// line via the associated messenger cmd /ZK/root/init before a run.
void rootStorageManager::CreateROOTObjects() {
  if (ROOTObjectsExist) {
    G4cout << "\nZK ANNOUNCEMENT : ROOT objects are presently are open and writeable! ZK cannot\n"
     <<   "                  initialize a new set until the existing ROOT objects are written to\n"
     <<   "                  to disk via the /ZK/root/write command.\n"
     << G4endl;

    return;
  }

  GenerateFileNames();

  // Create a ROOT TFile that will store all ROOT objects
  if (ROOTFile) delete ROOTFile;
  ROOTFile = new TFile(ROOTFileName, "recreate");

  // An example TTree to hold, for example, event-level data. The user
  // can create his/her own branches to store data. Note that data
  // members corresponding to data should be created such that branch
  // memory addresses can be correctly assigned
  ROOTTree = new TTree("ROOTTree", "TTree to hold event-level data");
  ROOTTree->Branch("energyDeposition", &energyDeposition);
  ROOTTree->Branch("weight", &weight);
  ROOTTree->Branch("beamEnergy", &beamEnergy);
  ROOTTree->Branch("gammaEnergy", &gammaEnergy);

  ROOTTreeFoil = new TTree("ROOTTreeFoil", "TTree to hold event-level data from foil");
  ROOTTreeFoil->Branch("energy",  &energy);
  ROOTTreeFoil->Branch("energyInitial", &energyInitial);
  ROOTTreeFoil->Branch("angle",   &angle );
  ROOTTreeFoil->Branch("weightF", &weightF);
  ROOTTreeFoil->Branch("isNRF",   &isNRF);
  ROOTTreeFoil->Branch("creatorProcessIndex", &creatorProcessIndex);
  ROOTTreeFoil->Branch("resonanceSample", &resonanceSample);

  ROOTTreeFoilInc = new TTree("ROOTTreeFoilInc", "TTree to tally photons stepping into foil");
  ROOTTreeFoilInc->Branch("foilEnergyInc", &foilEnergyInc);
  ROOTTreeFoilInc->Branch("foilWeightInc", &foilWeightInc);

  ROOTTreePlateInc = new TTree("ROOTTreePlateInc", "TTree to tally photons stepping into plate");
  ROOTTreePlateInc->Branch("plateEnergyInc", &plateEnergyInc);
  ROOTTreePlateInc->Branch("plateWeightInc", &plateWeightInc);

  ROOTObjectsExist = true;
}


// Method to write all ROOT objects to disk after a G4 run has
// completed. Note that this method should be called from the G4 cmd
// line via the associated messenger cmd /ZK/root/close after the
// run(s). This function is automatically called upon G4 session
// termination if ROOT objects still need to be written to prevent the
// user from losing data.
void rootStorageManager::WriteROOTObjects(G4bool EmergencyWrite) {
  if (!ROOTObjectsExist) {
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects presently do not exist and are, therefore, unwritable!\n"
           <<   "                  They should be created via the /ZK/root/init command before\n"
     <<   "                  events are recorded and then written to disk.\n"
     << G4endl;

    return;
  }

  if (EmergencyWrite)
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects that presently exist are being written to disk in\n"
           <<   "                  emergency fashion to avoid losing critical data before the simulation\n"
     <<   "                  terminates. Please issue the /ZK/root/write command before exiting\n"
     <<   "                  ZK to avoid this message.\n"
     << G4endl;

  // Output ROOT objects to the TFile
  ROOTTree->Write();
  ROOTTreeFoil->Write();
  ROOTTreeFoilInc->Write();
  ROOTTreePlateInc->Write();
  //ROOTFile->Close();

  G4cout << "Wrote data in ROOT objects to " << ROOTFileName << G4endl;

  ROOTObjectsExist = false;
}


// An example of a method that could be called from anywhere in the
// code to save data in the ROOT TTree. The user would do the following:
// 0. Pass this method the data to save
// 1. Save the data to a class member whose address has been registered
//    with the associated ROOT TTree branch
// 2. Call the TTree::Fill() function
void rootStorageManager::FillTree(G4double eD, G4double w, G4double eB, G4double eG) {
  // Save data to TTree here
  energyDeposition = eD;
  weight = w;
  beamEnergy = eB;
  gammaEnergy = eG;

  ROOTFile->cd();
  ROOTTree->Fill();
}

void rootStorageManager::FillTreeFoil(G4double e, G4double ei, G4double a, G4double w, G4bool n, G4int ci, G4bool rs) {
  energy  = e;
  energyInitial = ei;
  angle   = a;
  weightF = w;
  isNRF   = n;
  creatorProcessIndex = ci;
  resonanceSample = rs;

  ROOTFile->cd();
  ROOTTreeFoil->Fill();
}

void rootStorageManager::FillTreeFoilInc(G4double e, G4double w) {
  foilEnergyInc = e;
  foilWeightInc = w;
  ROOTFile->cd();
  ROOTTreeFoilInc->Fill();
}

void rootStorageManager::FillTreePlateInc(G4double e, G4double w) {
  plateEnergyInc = e;
  plateWeightInc = w;
  ROOTFile->cd();
  ROOTTreePlateInc->Fill();
}


// Method provided to generate ROOT TFile names. This is principally
// useful for batch processing so that each job can open its own
// TFile with a unique name based on its random seed.
void rootStorageManager::GenerateFileNames() {
  if (ROOTObjectsExist) {
    G4cout << "\nZK ANNOUNCEMENT : The ROOT objects presently exist and have already been assigned a name!\n"
     <<   "                  To change file names, please finish event processing and write the ROOT\n"
     <<   "                  objects to disk via the /ZK/root/write command. Then you are free\n"
     <<   "                  to set a new file name for subsequent event storage.\n"
     << G4endl;
    return;
  }


  // If the user has not specified a nonzero random seed / file ID,
  // the default file name uses the date/time as a unique identifier.
  std::stringstream ss;

  if (ROOTFileName == "") {
    if (find > 0) {
      char b2[256];
      snprintf(b2, sizeof(b2), "ZKExp_outfile_%d.root", find);
      G4String fname(b2);
      ROOTFileName = fname;
    } else {
      time_t theTime;
      struct tm *timeInfo;

      time(&theTime);
      timeInfo = localtime(&theTime);

      const int N = 100;
      char buffer[N];
      strftime(buffer, N, "ZK_%d%b%y_%H.%M", timeInfo);
      ss << buffer << ".root";
      ROOTFileName = ss.str();
    }
  }
}
