// Sensitive detector class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu
// borrows heavily from examples/medical/GammaTherapy/include/PhantomSD.hh

#include "sensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ParticleTypes.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Integrator.hh"

#include "rootStorageManager.hh"
#include "eventAction.hh"
#include "eventInformation.hh"

sensitiveDetector::sensitiveDetector(const G4String& name)
  : G4VSensitiveDetector(name),
  fCounter(0), totalEdep(0.), hits(0), incidentParticles(0),
  weight(0.), beamEnergy(0.), gammaEnergy(0.), massOfLogical(0.)
{}

sensitiveDetector::~sensitiveDetector()
{;}

void sensitiveDetector::Initialize(G4HCofThisEvent*) {
  ++fCounter;
  ResetTrackEdep();
}

void sensitiveDetector::ResetTrackEdep() {
  trackEdep = 0.0;
}

void sensitiveDetector::AccumulateTrackEdep(G4double edep) {
  trackEdep += edep;
}

void sensitiveDetector::AccumulateTotalEdep(G4double edep) {
  totalEdep += edep;
}

G4double sensitiveDetector::GetTrackEdep() {
  return trackEdep;
}

G4double sensitiveDetector::GetTotalEdep() {
  return totalEdep;
}

void sensitiveDetector::IncrementHits() {
  ++hits;
}

G4double sensitiveDetector::GetHits() {
  return hits;
}

void sensitiveDetector::IncrementIncidentParticles() {
  ++incidentParticles;
}

G4int sensitiveDetector::GetIncidentParticles() {
  return incidentParticles;
}

void sensitiveDetector::SetMassOfLogical(G4double m) {
  massOfLogical = m;
  G4cout << "Mass corresponding to " << this->GetName() << " set as " << massOfLogical/kg << " kg." << G4endl;
}

G4double sensitiveDetector::GetMassOfLogical() {
  return massOfLogical;
}


G4bool sensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory*) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep > 0.0) {
    G4ThreeVector p1 = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector p2 = aStep->GetPostStepPoint()->GetPosition();
    G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();

    // accumulate deposited energy for dose calculation, and increment hit counter
    eventInformation* info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    weight = info->GetWeight();
    AccumulateTrackEdep(edep * weight); // this weighting only works for sums, not individual tracks
    IncrementHits();
  }

  return true;
}


void sensitiveDetector::EndOfEvent(G4HCofThisEvent*) {
  G4double track_edep = GetTrackEdep();
  if (track_edep > 0.0) {
    // tally hits, energy deposition
    AccumulateTotalEdep(track_edep);
    IncrementIncidentParticles();

    // get the info of the (initial) event from the eventInformation
    eventInformation* info = (eventInformation*)(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    weight      = info->GetWeight();
    beamEnergy  = info->GetBeamEnergy();
    gammaEnergy = info->GetGammaEnergy();
  }
}


void sensitiveDetector::clear()
{}

void sensitiveDetector::PrintAll()
{}

void sensitiveDetector::PrintSummary() {
  G4cout << G4endl;
  G4cout << GetName() << " sensitive detector summary (for one job *only*):" << G4endl;
  G4cout << "  hits:                          " << GetHits()              << G4endl;
  G4cout << "  incident particles:            " << GetIncidentParticles() << G4endl;
  G4cout << "  total energy deposition (MeV): " << GetTotalEdep()/MeV     << G4endl;
  G4cout << "  total dose (Gy):               " << GetTotalEdep()/GetMassOfLogical() / (joule/kg) << G4endl;
  G4cout << G4endl;
}

