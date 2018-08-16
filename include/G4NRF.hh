//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// *                                                                  *
// * This version (C) 2018 Jayson Vavrek under the MIT License.       *
// ********************************************************************
//
//
// $Id: G4NRF.hh,v 1.1.1.1 2007/01/13 19:44:02 jordan Exp $
// GEANT4 tag $Name:  $
//
//      ------------ G4NRF physics process --------
//
// G4NRF Class description
//
// This class manages NRF, inheriting from G4VDiscreteProcess
//
// Class description - end

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4NRF_h
#define G4NRF_h 1

#include <fstream>
using std::ofstream;

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OrderedTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4NRFNuclearLevelManager.hh"
#include "AngularCorrelation.hh"

#include "G4VProcess.hh"

#include "G4Integrator.hh"
#include "c2_function.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4NRF : public G4VDiscreteProcess {
 public:
  G4NRF(const G4String& processName = "NRF", G4bool Verbose_in = false,
        G4bool use_xsec_tables_in = true, G4bool use_xsec_integration_in = true,
        G4bool force_isotropic_in = false);

  ~G4NRF();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void PrintInfoDefinition();

  G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);

  G4VParticleChange *PostStepDoIt(const G4Track& track, const G4Step&  step);

  void SetVerbose() {Verbose = true;}

  void SetParam_x(G4double x);
  G4double GetParam_x() const;

  void SetParam_t(G4double t);
  G4double GetParam_t() const;

  G4double expIntegrand(G4double y) const;
  G4double PsiIntegral(G4double x, G4double t, G4int nMeshpoints = 300, G4double sigmaBound = 4.0);
  G4double InterpolateCrossSection(const G4NRFNuclearLevel* pLevel, G4double GammaEnergy);

  void print_to_standalone(ofstream& file);

 private:
  G4NRF & operator=(const G4NRF &right);
  G4NRF(const G4NRF&);
  G4NRFNuclearLevelManager* pNuclearLevelManager;
  G4double NRF_xsec_calc_gaus(G4double GammaEnergy, const G4NRFNuclearLevel* pLevel);
  G4double NRF_xsec_calc(G4double GammaEnergy, G4NRFNuclearLevel* pLevel,
    G4bool useTables, G4int nMeshpoints = 300, G4double sigmaBound = 4.0);
  void MakeCrossSectionTable(G4NRFNuclearLevel *pLevel, G4double Teff);

  void SetupMultipolarityInfo(const G4int nLevel,
      const G4double E_gamma,
      const G4int jgamma,
      const G4NRFNuclearLevel* pLevel,
      const G4NRFNuclearLevel* pLevel_next,
      G4double& J0, G4double& J, G4double& Jf,
      G4int& L1, G4int& L2,
      G4double& Delta1, G4double& Delta2);

  void AssignMultipoles(const G4NRFNuclearLevel* pLevel,
      const G4double E_gamma,
      const G4int jgamma,
      const G4double Ji, const G4double Pi,
      const G4double Jf, const G4double Pf,
      G4int& L, G4double& Delta);

  G4ThreeVector SampleCorrelation(const G4double Ji, const G4double J,  const G4double Jf,
          const G4int L1, const G4int L2,
          const G4double Delta1, const G4double Delta2);

  G4ThreeVector SampleIsotropic();

  G4int FindMin_L(const G4double Ji, const G4double Pi,
    const G4double Jf, const G4double Pf, char& transition);

  G4bool ForbiddenTransition(const G4NRFNuclearLevel* pLevel,
    const G4NRFNuclearLevel* pLevel_next);

  G4double MixingRatio_WeisskopfEstimate(char multipole, const G4int L,
           const G4double E_gamma,
           const G4double Pi,
           const G4double Pf);

  G4double Lamda_Weisskopf(char multipole, const G4int L,
         const G4double E_gamma);

  G4int A_excited;
  G4int Z_excited;
  const G4NRFNuclearLevel* pLevel_excited;

  Angular_Correlation* pAngular_Correlation;

  G4bool Verbose;
  const G4bool use_xsec_tables;
  const G4bool use_xsec_integration;
  const G4bool force_isotropic_ang_corr;

  G4double param_x;
  G4double param_t;

  // numerical integration
  G4Integrator<const G4NRF, G4double(G4NRF::*)(G4double) const> integrator;
};

inline G4bool G4NRF::IsApplicable(const G4ParticleDefinition& particle) {
  return (&particle == G4Gamma::Gamma());
}

#endif

