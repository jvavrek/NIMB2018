// Physics list class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu

#include "physicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4DecayPhysics.hh"

#include "G4NRFPhysics.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"
//#include "G4EmLowEPPhysics.hh"

#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StoppingPhysics.hh"

#include "G4HadronicProcessStore.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"


physicsList::physicsList(G4bool use_xsec_tables_in, G4bool use_xsec_integration_in, G4bool force_isotropic_in)
  : use_xsec_tables(use_xsec_tables_in),
    use_xsec_integration(use_xsec_integration_in),
    force_isotropic(force_isotropic_in) {
  defaultCutValue = 0.05*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForPositron = defaultCutValue;
  cutForProton = defaultCutValue;

  G4HadronicProcessStore::Instance()->SetVerbose(0);

  ConstructPhysics();
}


physicsList::~physicsList()
{;}

void physicsList::ConstructPhysics() {
  G4DecayPhysics *theDecayPhysics = new G4DecayPhysics();
  theDecayPhysics->ConstructParticle();
  RegisterPhysics(theDecayPhysics);

  // Add NRF to the physicsList
  RegisterPhysics(new G4NRFPhysics("NRF", use_xsec_tables, use_xsec_integration, force_isotropic));
  G4cout << "\nAdded NRF to the physicsList.\n" << G4endl;

  // Add the rest of the usual suspects
  RegisterPhysics(new G4EmStandardPhysics());
  //RegisterPhysics(new G4EmLowEPPhysics());
  RegisterPhysics(new G4EmExtraPhysics(verboseLevel));

  // Ion stopping-in-matter physics
  RegisterPhysics( new G4StoppingPhysics(verboseLevel) );

  // Neutron tracking cuts for optimized simulation
  G4NeutronTrackingCut *theNeutronTrackingCut = new G4NeutronTrackingCut(verboseLevel);
  theNeutronTrackingCut->SetTimeLimit(10*microsecond);
  theNeutronTrackingCut->SetKineticEnergyLimit(0.01*eV);
  RegisterPhysics( theNeutronTrackingCut );
}


void physicsList::ConstructParticle() {
  G4Gamma::Definition();
  G4Electron::Definition();
  G4Positron::Definition();
}


void physicsList::ConstructProcess() {
  G4VModularPhysicsList::ConstructProcess();
}


void physicsList::SetCuts() {
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  SetCutValue(cutForProton, "proton");
}
