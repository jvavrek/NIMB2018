#ifndef sensitiveDetector_h
#define sensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;

class sensitiveDetector : public G4VSensitiveDetector {
 public:
  sensitiveDetector(const G4String&);
  virtual ~sensitiveDetector();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void PrintAll();

  void     AccumulateTotalEdep(G4double edep);
  void     AccumulateTrackEdep(G4double edep);
  void     ResetTrackEdep();
  G4double GetTrackEdep();
  G4double GetTotalEdep();
  void     IncrementHits();
  G4double GetHits();
  void     IncrementIncidentParticles();
  G4int    GetIncidentParticles();
  void     SetMassOfLogical(G4double m);
  G4double GetMassOfLogical();
  void     PrintSummary();

 private:
  sensitiveDetector & operator=(const sensitiveDetector &right);
  sensitiveDetector(const sensitiveDetector&);

  G4int fCounter;
  G4double trackEdep;
  G4double totalEdep;
  G4int hits;
  G4int incidentParticles;

  G4double weight;
  G4double beamEnergy;
  G4double gammaEnergy;

  G4double massOfLogical;
};

#endif
