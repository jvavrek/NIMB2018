// User stepping action class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu

#include "steppingAction.hh"

#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicalVolumeStore.hh"
#include "eventInformation.hh"
#include "rootStorageManager.hh"
#include "G4Gamma.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


steppingAction::steppingAction()
{;}


steppingAction::~steppingAction()
{;}


void steppingAction::UserSteppingAction(const G4Step *currentStep) {
  G4Track *currentTrack = currentStep->GetTrack();

  // track only gammas
  if (currentTrack->GetParticleDefinition() == G4Gamma::Definition()) {
    G4StepPoint *prevPoint = currentStep->GetPreStepPoint();
    G4StepPoint *postPoint = currentStep->GetPostStepPoint();

    G4VPhysicalVolume *prevV = prevPoint->GetPhysicalVolume();
    G4VPhysicalVolume *postV = postPoint->GetPhysicalVolume();

    // tally the gamma energy
    G4double gammaEnergy = currentTrack->GetKineticEnergy();

    if (postV != NULL) {
      // get the eventInformation
      eventInformation* info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());

      // also place this cut here in the steppingAction, not just stackingAction
      // cuts out downscatters that don't create new tracks, giving a ~25% speedup
      if (gammaEnergy < 1.7*MeV) currentTrack->SetTrackStatus(fStopAndKill);

      info->SetGammaEnergy(gammaEnergy);
      G4double beamEnergy = info->GetBeamEnergy();
      G4double weightF = info->GetWeight();
      G4bool isNRF = false;
      G4bool resSample = info->GetResonanceSample();
      G4int creatorIndex = info->GetCreatorProcessIndex();

      // throw a warning if energy conservation appears to be broken for gammas
      // small deviations have been observed; attribute to a likely mismatch in database energies
      if (gammaEnergy > beamEnergy) {
        G4cerr << G4endl;
        G4cerr << "Warning: gammaEnergy " << gammaEnergy << " MeV > beamEnergy " << beamEnergy << " MeV!!" << G4endl;
        G4cerr << "  energy difference (keV) = " << (gammaEnergy-beamEnergy)/keV << G4endl;
        G4cerr << "  creator process = " << currentTrack->GetCreatorProcess()->GetProcessName() << G4endl;
        G4cerr << "  material = " << currentTrack->GetMaterial()->GetName() << G4endl;

        // if the discrepancy is > 1.0 keV, throw an actual error
        if (gammaEnergy-beamEnergy > 1.0*keV) {
          G4cerr << "  energy difference > 1.0 keV. Aborting..." << G4endl;
          exit(111);
        } else {
          G4cerr << "  energy difference < 1.0 keV likely due to NRF database mismatch. Carrying on..." << G4endl;
        }
      }

      // get all photons leaving the foil, regardless of whether they were created in the foil or not
      G4String VName_prev = prevV->GetLogicalVolume()->GetName();
      G4String VName_post = postV->GetLogicalVolume()->GetName();
      if (VName_prev == "world_L" && VName_post == "theFoil_L") {
        rootStorageManager::GetInstance()->FillTreeFoilInc(gammaEnergy, weightF);
      }
      if (VName_prev == "world_L" && VName_post == "slabV_L") {
        //rootStorageManager::GetInstance()->FillTreePlateInc(gammaEnergy, weightF);
      }

      if (currentTrack->GetCreatorProcess() != 0) {
        // one more sanity check
        G4LogicalVolume *Vmoth_prev = prevV->GetMotherLogical();
        if (Vmoth_prev != NULL) {
          // check if leaving the foil and into the world volume; would like a more flexible solution than hardcoded string names
          if (VName_prev == "theFoil_L" && VName_post == "world_L") {
            // check if the photon was created by NRF -- this neglects Compton scatters that still leave the foil
            G4String CPName = currentTrack->GetCreatorProcess()->GetProcessName();
            if (CPName == "NRF") isNRF = true;

            G4ThreeVector pdir = currentTrack->GetMomentumDirection();
            G4double angle = beamDir.angle(pdir);
            //G4cout << gammaEnergy << " " << angle << " " << weightF << " " << isNRF << G4endl;
            if (angle > 3.0*pi/4.0 && gammaEnergy > 1.7*MeV) {
              rootStorageManager::GetInstance()->FillTreeFoil(gammaEnergy, beamEnergy, angle, weightF, isNRF, creatorIndex, resSample);
            }
          }
        }
      }
    }
  }
}
