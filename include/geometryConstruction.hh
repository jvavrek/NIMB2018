#ifndef geometryConstruction_hh
#define geometryConstruction_hh 1

#include <stdlib.h>

#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4PhysicalVolumeStore.hh"
#include "sensitiveDetector.hh"


class geometryConstruction : public G4VUserDetectorConstruction {
 public:
  geometryConstruction(G4int);
  ~geometryConstruction();

  G4VPhysicalVolume *Construct();

 private:
  G4Material *Vacuum;
  G4PhysicalVolumeStore *PVStore;

  const G4int geom_mode; // geometry option flag

  sensitiveDetector *Pu_BS_SD;
  sensitiveDetector *HE_BS_SD;
  sensitiveDetector *WU_BS_SD;

  G4double x_foil; // thickness of the foil mother volume
};

#endif
