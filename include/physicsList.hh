#ifndef physicsList_hh
#define physicsList_hh 1

#include "G4VModularPhysicsList.hh"


class physicsList: public G4VModularPhysicsList {
 public:
  physicsList(G4bool, G4bool, G4bool);
  ~physicsList();

  void ConstructParticle();
  void ConstructPhysics();
  void ConstructProcess();
  void SetCuts();

 private:
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;

  G4bool use_xsec_tables;
  G4bool use_xsec_integration;
  G4bool force_isotropic;
};

#endif

