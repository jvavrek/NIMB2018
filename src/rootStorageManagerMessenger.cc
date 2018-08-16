#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include "rootStorageManager.hh"
#include "rootStorageManagerMessenger.hh"


rootStorageManagerMessenger::rootStorageManagerMessenger(rootStorageManager *RSM)
  : theRSManager(RSM) {
  ZKDirectory = new G4UIdirectory("/ZK/");
  ZKDirectory->SetGuidance("Custom ZK commands");

  rootDirectory = new G4UIdirectory("/ZK/root/");
  rootDirectory->SetGuidance("Settings for SDE I/O via persistent ROOT data storage");

  rootFileNameCmd = new G4UIcmdWithAString("/ZK/root/setFileName", this);
  rootFileNameCmd->SetGuidance("Sets the name of the ROOT file. If a name is not set manually, the");
  rootFileNameCmd->SetGuidance("default file name will be '<date>_<time>.root'.");
  rootFileNameCmd->SetParameterName("Choice", false);

  rootInitCmd = new G4UIcmdWithoutParameter("/ZK/root/init", this);
  rootInitCmd->SetGuidance("Create and initialize the required ROOT objects for persistent storage.");
  rootInitCmd->SetGuidance("Note that this cmd may only be used when existing ROOT objects have been");
  rootInitCmd->SetGuidance("written to disk.");
  rootInitCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  rootWriteCmd = new G4UIcmdWithoutParameter("/ZK/root/write", this);
  rootWriteCmd->SetGuidance("Writes the ROOT objects to disk. Note that this cmd may only be used if");
  rootWriteCmd->SetGuidance("the ROOT objects have already been created and initialized.");
  rootWriteCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}


rootStorageManagerMessenger::~rootStorageManagerMessenger() {
  delete rootWriteCmd;
  delete rootInitCmd;
  delete rootFileNameCmd;
  delete rootDirectory;
}


void rootStorageManagerMessenger::SetNewValue(G4UIcommand *cmd, G4String newValue) {
  if (cmd == rootFileNameCmd)
    theRSManager->SetFileName(newValue);

  if (cmd == rootInitCmd)
    theRSManager->CreateROOTObjects();

  if (cmd == rootWriteCmd)
    theRSManager->WriteROOTObjects();
}
