#ifndef PGA_hh
#define PGA_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4Event.hh"

//#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "TH1D.h"
#include "TRandom1.h"

class PGA : public G4VUserPrimaryGeneratorAction {
 public:
  PGA(G4int, G4int);
  ~PGA();

  G4ParticleGun *GetSource() {return TheSource;}

  void GeneratePrimaries(G4Event *anEvent);

  G4double SampleUResonances();
  G4double SampleAlResonance();

  inline G4ThreeVector GetBeamDir() {return beamDir;}

 private:
  G4ParticleGun *TheSource;

  TRandom1 Random;

  G4ThreeVector beamDir;
  const G4int beam_mode;

  TH1D *hBrems;
  TH1D *hSample;
  TH1D *hBinary;
};

#endif
