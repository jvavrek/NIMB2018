// Event information class
// Based off of TrackInformation (http://geant4.slac.stanford.edu/Tips/event/1.html)
// Jayson Vavrek, MIT, 2016
// jvavrek@mit.edu

#include "eventInformation.hh"
#include "G4ios.hh"

//G4Allocator<eventInformation> anEventInformationAllocator;

eventInformation::eventInformation() {
  weight = 0.;
  beamEnergy = 0.;
  gammaEnergy = 0.; // might be better to have an init list

  creatorProcessIndex = 0;
}

eventInformation::eventInformation(const G4Event* anEvent) {
  weight = 0.; //aTrack->GetWeight();
  beamEnergy = 0.;
  gammaEnergy = 0.;

  creatorProcessIndex = 0;

  resSample = false;
}

eventInformation::eventInformation(const eventInformation* anEventInfo) {
  weight = anEventInfo->GetWeight();
  beamEnergy = anEventInfo->GetBeamEnergy();
  gammaEnergy = anEventInfo->GetGammaEnergy();

  creatorProcessIndex = anEventInfo->GetCreatorProcessIndex();
}

eventInformation::~eventInformation() {;}

void eventInformation::SetWeight(G4double x) {
  weight = x;
}

void eventInformation::SetBeamEnergy(G4double x) {
  beamEnergy = x;
}

void eventInformation::SetGammaEnergy(G4double x) {
  gammaEnergy = x;
}

void eventInformation::SetCreatorProcessIndex(G4int index) {
  creatorProcessIndex = index;
}

void eventInformation::IncrementCreatorProcessIndex(G4int inc) {
  creatorProcessIndex += inc;
}

void eventInformation::SetResonanceSample(G4bool res) {
  resSample = res;
}

void eventInformation::Print() const {
  //G4cout << "Original track ID " << originalTrackID << " at " << originalPosition << G4endl;
}
