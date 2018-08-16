#ifndef eventInformation_h
#define eventInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4Allocator.hh"
#include "G4VUserEventInformation.hh"

class eventInformation : public G4VUserEventInformation {
 public:
  eventInformation();
  eventInformation(const G4Event*);
  eventInformation(const eventInformation*);
  virtual ~eventInformation();

  inline G4double GetWeight() const {return weight;}
  void SetWeight(G4double);

  inline G4double GetBeamEnergy() const {return beamEnergy;}
  void SetBeamEnergy(G4double);

  inline G4double GetGammaEnergy() const {return gammaEnergy;}
  void SetGammaEnergy(G4double);

  inline G4int GetCreatorProcessIndex() const {return creatorProcessIndex;}
  void SetCreatorProcessIndex(G4int);
  void IncrementCreatorProcessIndex(G4int);

  inline G4bool GetResonanceSample() const {return resSample;}
  void SetResonanceSample(G4bool);

  void Print() const;

 private:
  G4double weight;
  G4double beamEnergy;
  G4double gammaEnergy;

  G4int creatorProcessIndex;

  G4bool resSample;
};

#endif
