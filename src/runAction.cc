#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

#include "runAction.hh"

#include "rootStorageManager.hh"
#include "geometryConstruction.hh"
#include "sensitiveDetector.hh"
#include "PGA.hh"


runAction::runAction()
{;}


runAction::~runAction()
{;}


void runAction::BeginOfRunAction(const G4Run *)
{;}


void runAction::EndOfRunAction(const G4Run *currentRun) {
  // print the incident particles and energy deposition
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  sensitiveDetector *Pu_BS_SD = (sensitiveDetector*) SDMan->FindSensitiveDetector("Pu_BS_SD", false); // warning = false if not found
  sensitiveDetector *HE_BS_SD = (sensitiveDetector*) SDMan->FindSensitiveDetector("HE_BS_SD", false);
  sensitiveDetector *WU_BS_SD = (sensitiveDetector*) SDMan->FindSensitiveDetector("WU_BS_SD", false);

  if (WU_BS_SD) WU_BS_SD->PrintSummary();
  if (HE_BS_SD) HE_BS_SD->PrintSummary();
  if (Pu_BS_SD) Pu_BS_SD->PrintSummary();

  if (WU_BS_SD || HE_BS_SD || Pu_BS_SD) {
    G4cout << "Also keep in mind that the dose will be affected by the following:" << G4endl;
    G4cout << "  - any cuts in the stepping and stackingAction classes" << G4endl;
    G4cout << "  - the total number of particles run" << G4endl;
    G4cout << "  - whether the sampling cuts off low-energy photons" << G4endl;
    G4cout << "  - whether or not the weighting of Edep was done properly" << G4endl;
    G4cout << G4endl;
  }
}
