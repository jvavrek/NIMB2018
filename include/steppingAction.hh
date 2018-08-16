#ifndef steppingAction_hh
#define steppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

class steppingAction : public G4UserSteppingAction {
 public:
  steppingAction();
  ~steppingAction();

  void UserSteppingAction(const G4Step *);

  inline void SetBeamDir(G4ThreeVector bd) {beamDir = bd;}

 private:
  G4ThreeVector beamDir;
};

#endif
