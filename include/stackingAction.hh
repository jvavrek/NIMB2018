#ifndef stackingAction_hh
#define stackingAction_hh 1

#include "G4UserStackingAction.hh"
#include "G4Track.hh"

class stackingAction : public G4UserStackingAction {
 public:
  stackingAction();
  ~stackingAction();

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

  inline void SetBeamDir(G4ThreeVector bd) {beamDir = bd;}

 private:
  G4ThreeVector beamDir;
};

#endif
