// User primary generator action class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu

#include "PGA.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4Gamma.hh"

#include "TFile.h"
#include "TROOT.h"

#include "eventInformation.hh"


// Create the bremsstrahlung or nearly-monochromatic pencil beam source.
PGA::PGA(G4int randomSeed, G4int beam_mode_in)
  : beam_mode(beam_mode_in) {
  // validate the random seed mode
  if (randomSeed < 0) {
    G4cout << "Error! Invalid randomSeed of " << randomSeed << " in PGA." << G4endl;
    G4cout << "Aborting..." << G4endl;
    exit(46);
  }

  // validate the beam mode
  if (beam_mode < 0 || beam_mode > 2) {
    G4cout << "Error! Invalid beam_mode of " << beam_mode << " in PGA." << G4endl;
    G4cout << "Aborting..." << G4endl;
    exit(47);
  }

  gRandom->SetSeed(randomSeed);

  G4cout << "PGA settings:" << G4endl;
  G4cout << "  randomSeed = " << randomSeed << G4endl;
  G4cout << "  beam_mode  = " << beam_mode  << G4endl;

  // establish the basics of the source
  TheSource = new G4ParticleGun();
  TheSource->SetParticleDefinition(G4Gamma::Definition());

  beamDir = G4ThreeVector(1, 0, 0); // downbeam direction is along x, or (1, 0, 0)

  if (beam_mode == 2) {
    // file contains the normalized brems distribution p(E), sampling distribution s(E),
    // and binary 0/1 for off/on resonance useful in weighting
    TFile *fin = TFile::Open("brems_distributions.root");
    hBrems  = (TH1D*) fin->Get("hBrems");
    hSample = (TH1D*) fin->Get("hSample");
    hBinary = (TH1D*) fin->Get("hBinary");

    if (hBrems && hSample && hBinary) {
      G4cout << "Imported brems and sampling distributions from " << fin->GetName() << G4endl << G4endl;
    } else {
      G4cout << "Error reading from file " << fin->GetName() << G4endl;
      exit(1);
    }
  }
}


PGA::~PGA() {
  delete TheSource;
}


void PGA::GeneratePrimaries(G4Event *currentEvent) {
  G4double eSample;
  //eSample = Random.Uniform(2.0*MeV, 2.5*MeV);
  //eSample = 2.2*MeV;
  //eSample = 2.20901100546*MeV;
  //eSample = 2.17601067911*MeV;
  //eSample = 2.21210728333*MeV;
  //eSample = 1.73354686425*MeV;
  //eSample = 2.00336916735*MeV;
  //eSample = 1.81525753275*MeV;
  //eSample = 2.423*MeV;
  //eSample = 2.43171328057*MeV;
  //eSample = 2.43321324158*MeV-2*eV;
  //eSample = 2.57751485869*MeV;
  //eSample = 2.56641473101*MeV;
  //eSample = 1.78200716196*MeV;
  //eSample = 1.84600768563*MeV;
  //eSample = 2.1760106791*MeV;
  //eSample = 2.21210728335*MeV;

  if (beam_mode == 0) {
    eSample = SampleUResonances(); // sample the U-238 resonances near 2.176, 2.209, 2.245 MeV
  } else if (beam_mode == 1) {
    eSample = SampleAlResonance(); // sample the primary Al-27 resonance near 2.212 MeV
  } else if (beam_mode == 2) {
    eSample = hSample->GetRandom()*MeV; // sample the resonances specified by hSample
  } else {
    G4cerr << "Error! Invalid beam_mode = " << beam_mode << G4endl;
    G4cerr << "Aborting..." << G4endl;
    exit(40);
  }

  G4double ySample = 0.0;
  G4double zSample = 0.0;

  TheSource->SetParticleEnergy(eSample);
  TheSource->SetParticlePosition(G4ThreeVector(-1.0*m, ySample, zSample));
  TheSource->SetParticleMomentumDirection(beamDir);

  TheSource->GeneratePrimaryVertex(currentEvent);

  // Calculate the weights. For the simple simulations, weights are 1. For the brems
  // simulations, weights are from importance sampling.
  G4double w;
  if (beam_mode < 2) {
    w = 1.0;
  } else {
    G4double    s = hSample->GetBinContent(hSample->GetXaxis()->FindBin(eSample));
    G4double dNdE = hBrems->GetBinContent(hBrems->GetXaxis()->FindBin(eSample));
    w = dNdE/s;
  }

  // if histogram exists, find whether it was a resonance sample or not
  G4bool onRes;
  if (hBinary)
    onRes = hBinary->GetBinContent(hBinary->GetXaxis()->FindBin(eSample));
  else
    onRes = false;

  // pass the event information up the chain
  eventInformation* anInfo = new eventInformation(currentEvent);
  anInfo->SetWeight(w);
  anInfo->SetBeamEnergy(eSample);
  anInfo->SetResonanceSample(onRes);
  currentEvent->SetUserInformation(anInfo);
}


G4double PGA::SampleUResonances() {
  std::vector<double> er;
  er.push_back(2.176*MeV);
  er.push_back(2.209*MeV);
  er.push_back(2.245*MeV);

  G4int idx = Random.Integer(er.size());
  G4double de = 25.0*eV;

  return Random.Uniform(er[idx]-de, er[idx]+de);
}

G4double PGA::SampleAlResonance() {
  G4double de = 25.0*eV;
  G4double er = 2.212107*MeV;
  return Random.Uniform(er-de, er+de);
}

