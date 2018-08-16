// User stacking action class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu

#include "stackingAction.hh"

#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "rootStorageManager.hh"
#include "eventInformation.hh"
#include "PGA.hh"
#include "G4RunManager.hh"


stackingAction::stackingAction()
{;}


stackingAction::~stackingAction()
{;}


G4ClassificationOfNewTrack stackingAction::ClassifyNewTrack(const G4Track* currentTrack) {
  // apply cuts in angle while roughly behind the detector
  G4ThreeVector pdir = currentTrack->GetMomentumDirection();
  G4double angle = beamDir.angle(pdir);
  G4double trackX = currentTrack->GetPosition().x();
  G4double trackE = currentTrack->GetTotalEnergy();

  // these might be useful later:
  // get the user eventInformation and the track creator process
  //const G4VProcess *cproc = currentTrack->GetCreatorProcess();
  eventInformation* info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
  G4ParticleDefinition *pdef = currentTrack->GetDefinition();

  // if a new track is created in the object with large angle and/or low energy, stop tracking it
  if ((angle > pi/4.0 || trackE < 1.0*MeV) && trackX < 23.5*cm) return fKill;

  // kill neutrons (probably not important)
  if (pdef == G4Neutron::Definition()) return fKill;

  // kill all photons below 1.7 MeV. Also in steppingAction
  if (pdef == G4Gamma::Definition() && trackE < 1.7*MeV) return fKill;

  // kill any secondaries (e.g. Compton/photo/pair e-) in the object
  G4int trackID = currentTrack->GetTrackID();
  if (trackID != 1 && trackX < 23.5*cm) return fKill;

  // tag continuum photons with their production method
  // change bool below to false to forgo tagging if unnecessary or too slow
  const G4bool continuumTagging = true;
  if (continuumTagging) {
    // tag via the secondary electrons that give off the brem
    // probably not very robust, since one e+/e- track sets the index for the whole event
    if (trackID != 1 && trackX > 23.5*cm && (pdef == G4Electron::Definition() || pdef == G4Positron::Definition())) {
      G4String cpname = currentTrack->GetCreatorProcess()->GetProcessName();
      if      (cpname == "phot" ) info->IncrementCreatorProcessIndex( 1); // photoelectric absorption -- e- brem
      else if (cpname == "compt") info->IncrementCreatorProcessIndex(10); // Compton scattering       -- e- brem
      else if (cpname == "conv" ) info->IncrementCreatorProcessIndex(50); // pair production          -- e+/- brem plus e+ inflight annihil
    }
  }

  return fUrgent;
}
