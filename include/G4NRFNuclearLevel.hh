//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------------
//
//      File name:     G4NRFNuclearLevel
//      Author:        David Jordan (david.jordan@pnl.gov)
//      Creation date: October 2006
//
//      Adapted from:  G4NuclearLevel,
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//      Creation date: 24 October 1998
//
//      Modifications:
//        14 Oct. 2006, David Jordan (DavidJordan@pnl.gov)
//              Adapted G4NuclearLevel class to NRF physics process.
//
//        G4NuclearLevel history follows:
//
//        09 Sep. 2002, Fan Lei  (flei@space.qinetiq.com)
//              Added IC probability when calculate the channel probabilities in
//              MakeProbabilities().
//
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions.
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data.
//
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVEL_HH
#define G4NUCLEARLEVEL_HH

#include <vector>
#include <fstream>

#include "globals.hh"
#include "c2_function.hh"

using std::ofstream;

class G4NRFNuclearLevelManager;

class G4NRFNuclearLevel {
 public:
  G4NRFNuclearLevel(const G4int nLevel, const G4int nucleusZ, const G4int nucleusA,
     const G4double E_1stExcitedState,
     const G4double energy, const G4double halfLife,
     const G4double angularMomentum, const std::vector<double>& eGamma,
     const std::vector<double>& wGamma, const std::vector<double>& polarities,
     const std::vector<double>& kCC, const std::vector<double>& l1CC,
     const std::vector<double>& l2CC, const std::vector<double>& l3CC,
     const std::vector<double>& m1CC, const std::vector<double>& m2CC,
     const std::vector<double>& m3CC, const std::vector<double>& m4CC,
     const std::vector<double>& m5CC, const std::vector<double>& nPlusCC,
     const std::vector<double>& totalCC, const G4bool Verbose = false
     );

  ~G4NRFNuclearLevel();

  const std::vector<double>& GammaEnergies() const;

  const std::vector<double>& GammaWeights() const;

  const std::vector<double>& GammaProbabilities() const;

  const std::vector<double>& GammaCumulativeProbabilities() const;

  const std::vector<double>& GammaPolarities() const;

  const std::vector<double>& KConvertionProbabilities() const;

  const std::vector<double>& L1ConvertionProbabilities() const;

  const std::vector<double>& L2ConvertionProbabilities() const;

  const std::vector<double>& L3ConvertionProbabilities() const;

  const std::vector<double>& M1ConvertionProbabilities() const;

  const std::vector<double>& M2ConvertionProbabilities() const;

  const std::vector<double>& M3ConvertionProbabilities() const;

  const std::vector<double>& M4ConvertionProbabilities() const;

  const std::vector<double>& M5ConvertionProbabilities() const;

  const std::vector<double>& NPlusConvertionProbabilities() const;

  const std::vector<double>& TotalConvertionProbabilities() const;

  const std::vector<int>&    MultipoleNumModes() const;

  const std::vector<char>&   MultipoleMode1() const;

  const std::vector<int>&    MultipoleL1() const;

  const std::vector<char>&   MultipoleMode2() const;

  const std::vector<int>&    MultipoleL2() const;

  const std::vector<double>& MultipoleMixingRatio() const;

  const std::vector<int>&    MultipoleMixRatioSignFlag() const;

  G4int nLevel() const;

  inline G4double Energy() const {return _energy;}

  G4int Z() const;

  G4int A() const;

  G4double AngularMomentum() const;

  G4double Parity() const;

  G4double HalfLife() const;

  G4int NumberOfGammas() const;

  G4double Width() const;

  G4double Width0() const;

  G4double SelectGamma(G4int& igamma) const;

  G4double MaxGammaEnergy() const;

  void PrintAll() const;

  void PrintAllTabular(ofstream& file, double gsSpin, double TDebye) const;

  G4bool GetInvalidLevel();
  void SetInvalidLevel();

  const interpolating_function_p<G4double>* GetCrossSectionTable() const;
  void SetCrossSectionTable(interpolating_function_p<G4double> *f);

  G4bool operator==(const G4NRFNuclearLevel &right) const;
  G4bool operator!=(const G4NRFNuclearLevel &right) const;
  G4bool operator<(const G4NRFNuclearLevel &right) const;

  const G4NRFNuclearLevel& operator=(const G4NRFNuclearLevel &right) {
    if (this != &right) {
      _energies   = right._energies;
      _weights    = right._weights;
      _prob       = right._prob;
      _cumProb    = right._cumProb;
      _polarities = right._polarities;
      _kCC        = right._kCC;
      _l1CC       = right._l1CC;
      _l2CC       = right._l2CC;
      _l3CC       = right._l3CC;
      _m1CC       = right._m1CC;
      _m2CC       = right._m2CC;
      _m3CC       = right._m3CC;
      _m4CC       = right._m4CC;
      _m5CC       = right._m5CC;
      _nPlusCC    = right._nPlusCC;
      _totalCC    = right._totalCC;
      _Num_multipole              = right._Num_multipole;
      _Multipole_mode1            = right._Multipole_mode1;
      _Multipole_L1               = right._Multipole_L1;
      _Multipole_mode2            = right._Multipole_mode2;
      _Multipole_L2               = right._Multipole_L2;
      _Multipole_mixing_ratio     = right._Multipole_mixing_ratio;
      _Multipole_mixing_sign_flag = right._Multipole_mixing_sign_flag;

      _nLevel             = right._nLevel;
      _nucleusZ           = right._nucleusZ;
      _nucleusA           = right._nucleusA;
      _energy             = right._energy;
      _halfLife           = right._halfLife;
      _angularMomentum    = right._angularMomentum;
      _parity             = right._parity;
      _nGammas            = right._nGammas;
      _Tau                = right._Tau;
      _Tau0               = right._Tau0;
      _E_1stExcitedState  = right._E_1stExcitedState;
      _Ewidth_gamma       = right._Ewidth_gamma;
      _Ewidth_gamma0      = right._Ewidth_gamma0;
      _Ewidth_p           = right._Ewidth_p;
      _Ewidth_n           = right._Ewidth_n;
      _Ewidth_alpha       = right._Ewidth_alpha;
      _Verbose            = right._Verbose;
    }

    return *this;
  }

  G4NRFNuclearLevel(const G4NRFNuclearLevel &right) {
    if (this != &right) *this = right;
  }

 private:
  G4NRFNuclearLevel() {G4cout << "Calling default constructor" << G4endl;}

  void MakeProbabilities();
  void MakeCumProb();
  void RefreshWidth();
  void RefreshGammas();
  void MakeWidth0();
  void Calc_Tau0();

  G4bool Identify_GS_Transition();

  G4double GetGSProb() const;

  G4int Increment(G4int aF);

  std::vector<double> _energies;
  std::vector<double> _weights;
  std::vector<double> _prob;
  std::vector<double> _cumProb;
  std::vector<double> _polarities;
  std::vector<double> _kCC;
  std::vector<double> _l1CC;
  std::vector<double> _l2CC;
  std::vector<double> _l3CC;
  std::vector<double> _m1CC;
  std::vector<double> _m2CC;
  std::vector<double> _m3CC;
  std::vector<double> _m4CC;
  std::vector<double> _m5CC;
  std::vector<double> _nPlusCC;
  std::vector<double> _totalCC;
  std::vector<int>    _Num_multipole;
  std::vector<char>   _Multipole_mode1;
  std::vector<int>    _Multipole_L1;
  std::vector<char>   _Multipole_mode2;
  std::vector<int>    _Multipole_L2;
  std::vector<double> _Multipole_mixing_ratio;
  std::vector<int>    _Multipole_mixing_sign_flag;

  G4int    _nLevel;
  G4int    _nucleusZ;
  G4int    _nucleusA;
  G4double _energy;
  G4double _halfLife;
  G4double _angularMomentum;
  G4double _parity;
  G4int    _nGammas;
  G4double _Tau;
  G4double _Tau0;

  G4double _E_1stExcitedState;
  G4double _Ewidth_gamma;
  G4double _Ewidth_gamma0;
  G4double _Ewidth_p;
  G4double _Ewidth_n;
  G4double _Ewidth_alpha;

  G4bool   _Verbose;
  G4bool   invalidLevel;

  // cross section table for interpolation
  interpolating_function_p<G4double> *_cross_sec_interp_func;
};

#endif
