// Geometry construction class for the ZKExp NRF simulations
// Jayson Vavrek, MIT, 2015
// jvavrek@mit.edu

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4PVParameterised.hh"

#include "geometryConstruction.hh"
#include "rootStorageManager.hh"

geometryConstruction::geometryConstruction(G4int geom_mode_in)
  : geom_mode(geom_mode_in) {
  Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  // validate the geometry mode
  if (geom_mode < 0 || geom_mode > 3) {
    G4cout << "Error! Invalid geom_mode of " << geom_mode << " in geometryConstruction." << G4endl;
    G4cout << "Aborting..." << G4endl;
    exit(45);
  }

  G4cout << "geometryConstruction settings:" << G4endl;
  G4cout << "  geom_mode = " << geom_mode << G4endl;
}


geometryConstruction::~geometryConstruction()
{;}


G4VPhysicalVolume *geometryConstruction::Construct() {
  //-------------------------------------------------------------------------
  // Colours
  //                           RED GREEN BLUE
  const G4Colour  RED         (1.0, 0.0, 0.0);
  const G4Colour  GREEN       (0.0, 1.0, 0.0);
  const G4Colour  BLUE        (0.0, 0.0, 1.0);
  const G4Colour  YELLOW      (1.0, 1.0, 0.0);
  const G4Colour  CYAN        (0.0, 1.0, 1.0);
  const G4Colour  PURPLE      (1.0, 0.0, 1.0);
  const G4Colour  ORANGE      (1.0, 0.5, 0.0);
  const G4Colour  GRAY        (0.5, 0.5, 0.5);
  const G4Colour  LIGHTGRAY   (0.8, 0.8, 0.8);
  const G4Colour  DARKGRAY    (0.3, 0.3, 0.3);
  const G4Colour  DARKERGRAY  (0.1, 0.1, 0.1);
  const G4Colour  DARKRED     (0.5, 0.0, 0.0);
  const G4Colour  DARKGREEN   (0.0, 0.5, 0.0);
  const G4Colour  DARKBLUE    (0.0, 0.0, 0.5);
  const G4Colour  DARKYELLOW  (0.5, 0.5, 0.0);
  const G4Colour  DARKCYAN    (0.0, 0.5, 0.5);
  const G4Colour  DARKPURPLE  (0.5, 0.0, 0.5);
  const G4Colour  DARKORANGE  (0.5, 0.2, 0.0);
  const G4Colour  BLACK       (0.0, 0.0, 0.0);
  const G4Colour  WHITE       (1.0, 1.0, 1.0);
  const G4Colour  LIGHTBLUE   (0.5, 0.8, 1.0);

  //-------------------------------------------------------------------------
  // Materials

  G4double a;            // atomic mass
  G4double z;            // atomic number
  G4int iz;              // atomic number as integer
  G4int n;               // number of neutrons
  G4double density;      // density of material
  G4String name;         // name of element, isotope or material
  G4String symbol;       // symbol of element
  G4int ncomponents;     // number of components in element
  G4int natoms;          // number of atoms in weighting of material
  G4double abundance;    // fractional abundance
  G4double fractionmass; // fraction of material mass

  // Natural aluminum
  G4Isotope *Al27  = new G4Isotope(name = "Al27", iz = 13, n = 27, a = 26.982*g/mole);
  G4Element *NatAl = new G4Element(name = "Natural Aluminum", symbol = "NatAl", ncomponents = 1);
  NatAl->AddIsotope(Al27, abundance = 100*perCent);
  G4Material *NatAl_mat = new G4Material(name = "Natural Aluminum", density = 2.6989*g/cm3, ncomponents = 1);
  NatAl_mat->AddElement(NatAl, fractionmass = 100.0*perCent);

  // Natural lead
  G4Isotope *Pb204 = new G4Isotope(name = "Pb204", iz = 82, n = 204, a = 203.97304*g/mole);
  G4Isotope *Pb206 = new G4Isotope(name = "Pb206", iz = 82, n = 206, a = 205.97447*g/mole);
  G4Isotope *Pb207 = new G4Isotope(name = "Pb207", iz = 82, n = 207, a = 206.97590*g/mole);
  G4Isotope *Pb208 = new G4Isotope(name = "Pb208", iz = 82, n = 208, a = 207.97665*g/mole);
  G4Element *NatPb = new G4Element(name = "Natural Lead", symbol = "Pb", ncomponents = 4);
  NatPb->AddIsotope(Pb204, abundance =  1.5*perCent);
  NatPb->AddIsotope(Pb206, abundance = 23.6*perCent);
  NatPb->AddIsotope(Pb207, abundance = 22.6*perCent);
  NatPb->AddIsotope(Pb208, abundance = 52.3*perCent);

  G4Material* NatPb_mat = new G4Material(name = "Natural Lead", density = 11.35*g/cm3, ncomponents = 1);
  NatPb_mat->AddElement(NatPb, fractionmass = 100.0*perCent);

  // Natural iron
  G4Isotope *Fe54 = new G4Isotope(name = "Fe54", iz = 26, n = 54, a = 53.940*g/mole);
  G4Isotope *Fe56 = new G4Isotope(name = "Fe56", iz = 26, n = 56, a = 55.935*g/mole);
  G4Isotope *Fe57 = new G4Isotope(name = "Fe57", iz = 26, n = 57, a = 56.935*g/mole);
  G4Isotope *Fe58 = new G4Isotope(name = "Fe58", iz = 26, n = 58, a = 57.933*g/mole);
  G4Element *NatFe = new G4Element(name = "Natural Iron", symbol = "Fe", ncomponents = 4);
  NatFe->AddIsotope(Fe54, abundance =  5.8*perCent);
  NatFe->AddIsotope(Fe56, abundance = 91.8*perCent);
  NatFe->AddIsotope(Fe57, abundance =  2.1*perCent);
  NatFe->AddIsotope(Fe58, abundance =  0.3*perCent);

  G4Material* NatFe_mat = new G4Material(name = "Natural Iron", density = 7.874*g/cm3, ncomponents = 1);
  NatFe_mat->AddElement(NatFe, fractionmass = 100.0*perCent);

  // Natural hydrogen
  G4Isotope *H1 = new G4Isotope(name = "H1", iz = 1, n = 1, a = 1.01*g/mole);
  G4Isotope *H2 = new G4Isotope(name = "H2", iz = 1, n = 2, a = 2.01*g/mole);
  G4Element *NatH = new G4Element(name = "Natural Hydrogen", symbol = "NatH", ncomponents = 2);
  NatH->AddIsotope(H1, abundance = 0.9999);
  NatH->AddIsotope(H2, abundance = 0.0001);

  // Natural carbon
  G4Isotope *C12 = new G4Isotope(name = "C12", iz = 6, n = 12, a = 12.000*g/mole);
  G4Isotope *C13 = new G4Isotope(name = "C13", iz = 6, n = 13, a = 13.003*g/mole);
  G4Element *NatC = new G4Element(name = "Natural Carbon", symbol = "C", ncomponents = 2);
  NatC->AddIsotope(C12, abundance = 98.89*perCent);
  NatC->AddIsotope(C13, abundance =  1.11*perCent);

  // Natural nitrogen
  G4Isotope *N14 = new G4Isotope(name = "N14", iz = 7, n = 14, a = 14.003*g/mole);
  G4Isotope *N15 = new G4Isotope(name = "N15", iz = 7, n = 15, a = 15.000*g/mole);
  G4Element *NatN = new G4Element(name = "Natural Nitrogen", symbol = "N", ncomponents = 2);
  NatN->AddIsotope(N14, abundance = 99.64*perCent);
  NatN->AddIsotope(N15, abundance =  0.36*perCent);

  // Natural Oxygen
  G4Isotope *O16 = new G4Isotope(name = "O16", iz = 8, n = 16, a = 15.995*g/mole);
  G4Isotope *O17 = new G4Isotope(name = "O17", iz = 8, n = 17, a = 16.999*g/mole);
  G4Isotope *O18 = new G4Isotope(name = "O18", iz = 8, n = 18, a = 17.992*g/mole);
  G4Element *NatO = new G4Element(name = "Natural Oxygen", symbol = "O", ncomponents = 3);
  NatO->AddIsotope(O16, abundance = 99.757*perCent);
  NatO->AddIsotope(O17, abundance =  0.038*perCent);
  NatO->AddIsotope(O18, abundance =  0.205*perCent);

  // Uranium isotopes
  G4Isotope *U235 = new G4Isotope(name = "U235", iz = 92, n = 235, a = 235.044*g/mole);
  G4Isotope *U238 = new G4Isotope(name = "U238", iz = 92, n = 238, a = 238.051*g/mole);

  // Simple uranium-238
  G4Element *PureU238 = new G4Element(name = "PureU238", symbol = "PureU238", ncomponents = 1);
  PureU238->AddIsotope(U238, abundance = 1.00);
  G4Material *U238_simple = new G4Material(name = "Simple U-238", density = 19.052*g/cm3, ncomponents = 1);
  U238_simple->AddElement(PureU238, fractionmass = 100.0*perCent);

  // Simple uranium-235
  G4Element *PureU235 = new G4Element(name = "PureU235", symbol = "PureU235", ncomponents = 1);
  PureU235->AddIsotope(U235, abundance = 1.00);
  G4Material *U235_simple = new G4Material(name = "Simple U-235", density = 18.811*g/cm3, ncomponents = 1);
  U235_simple->AddElement(PureU235, fractionmass = 100.0*perCent);

  // Simple blend of 235 and 238, 95 (wt.%) enriched
  G4Material *Uenriched95 = new G4Material(name = "95pcw enriched uranium", density = 18.7*g/cm3, ncomponents = 2);
  Uenriched95->AddElement(PureU238, fractionmass =  5.0*perCent);
  Uenriched95->AddElement(PureU235, fractionmass = 95.0*perCent);

  // Plutonium isotopes
  G4Isotope *Pu239 = new G4Isotope(name = "Pu239", iz = 94, n = 239, a = 239.052*g/mole);
  G4Isotope *Pu240 = new G4Isotope(name = "Pu240", iz = 94, n = 240, a = 240.054*g/mole);

  // Simple plutonium-239
  G4Element *PurePu239 = new G4Element(name = "PurePu239", symbol = "PurePu239", ncomponents = 1);
  PurePu239->AddIsotope(Pu239, abundance = 1.00);
  G4Material *Pu239_simple = new G4Material(name = "Simple Pu-239", density = 19.41*g/cm3 , ncomponents = 1);
  Pu239_simple->AddElement(PurePu239, fractionmass = 100.0*perCent);

  // Simple plutonium-240
  G4Element *PurePu240 = new G4Element(name = "PurePu240", symbol = "PurePu240", ncomponents = 1);
  PurePu240->AddIsotope(Pu240, abundance = 1.00);
  G4Material *Pu240_simple = new G4Material(name = "Simple Pu-240", density = 19.496*g/cm3 , ncomponents = 1);
  Pu240_simple->AddElement(PurePu240, fractionmass = 100.0*perCent);

  // Simple definitions of WgPu and FgPu
  G4Material *WgPu_simple = new G4Material(name = "Simple WgPu", density = 19.84*g/cm3, ncomponents = 2);
  WgPu_simple->AddElement(PurePu239, fractionmass = 0.94);
  WgPu_simple->AddElement(PurePu240, fractionmass = 0.06);

  G4Material *FgPu_simple = new G4Material(name = "Simple FgPu", density = 19.84*g/cm3, ncomponents = 2);
  FgPu_simple->AddElement(PurePu239, fractionmass = 0.86);
  FgPu_simple->AddElement(PurePu240, fractionmass = 0.14);

  // HMX (high explosive)
  G4Material* HMX = new G4Material(name = "HMX", density = 1.890*g/cm3, ncomponents = 4);
  HMX->AddElement(NatH, fractionmass = 0.027227);
  HMX->AddElement(NatC, fractionmass = 0.162222);
  HMX->AddElement(NatN, fractionmass = 0.378361);
  HMX->AddElement(NatO, fractionmass = 0.432190);

  // composite foil of 25% by mass each of U-235, U-238, Pu-239, Pu-240; arbitrarily 19 g/cm3
  G4Material *composite25pc = new G4Material(name = "composite25pc", density = 19*g/cm3, ncomponents = 4);
  composite25pc->AddElement(PureU235,  fractionmass = 0.25);
  composite25pc->AddElement(PureU238,  fractionmass = 0.25);
  composite25pc->AddElement(PurePu239, fractionmass = 0.25);
  composite25pc->AddElement(PurePu240, fractionmass = 0.25);


  //-------------------------------------------------------------------------
  // Sensitive detector(s) initialization

  G4SDManager* SDMan = G4SDManager::GetSDMpointer();

  // three SDs for calculation of dose to the measurement object
  if (geom_mode == 2) {
    Pu_BS_SD = new sensitiveDetector("Pu_BS_SD");
    HE_BS_SD = new sensitiveDetector("HE_BS_SD");
    WU_BS_SD = new sensitiveDetector("WU_BS_SD");
  }


  //-------------------------------------------------------------------------
  // Volumes

  // Physical volume store
  PVStore = G4PhysicalVolumeStore::GetInstance();

  //----World volume
  //
  G4double worldX = 2*m;
  G4double worldY = 2*m;
  G4double worldZ = 2*m;

  G4Box *world_S = new G4Box("world_S", worldX, worldY, worldZ);

  G4LogicalVolume *world_L = new G4LogicalVolume(world_S,
             Vacuum,
             "world_L");

  G4VPhysicalVolume *world_P = new G4PVPlacement(0,
             G4ThreeVector(),
             world_L,
             "TheWorld",
             0,
             false,
             0);

  G4VisAttributes *worldVisAtt = new G4VisAttributes();
  worldVisAtt->SetVisibility(false);
  world_L->SetVisAttributes(worldVisAtt);

  // Shared attributes sets across the geometries
  G4VisAttributes *PuVisAtt = new G4VisAttributes(true, DARKORANGE);
  G4VisAttributes * UVisAtt = new G4VisAttributes(true, DARKBLUE);
  G4VisAttributes *AlVisAtt = new G4VisAttributes(true, LIGHTGRAY);
  G4VisAttributes *HEVisAtt = new G4VisAttributes(true, DARKGREEN);


  //----encrypting foil, homogenized mixture of isotopes
  //
  x_foil = 2.0*cm; // foil thickness (X in the glass equation)
  G4double y_foil = 50.0*cm; // foil height
  G4double z_foil = 50.0*cm; // foil width
  G4double d_foil = 1.00*m;  // distance from origin (centre of object) to centre of foil

  G4Box *theFoil_S = new G4Box("theFoil_S", x_foil/2.0, y_foil/2.0, z_foil/2.0);

  G4Material *theFoil_mat = 0;
  if (geom_mode == 0)
    theFoil_mat = U238_simple;
  else if (geom_mode == 1)
    theFoil_mat = NatAl_mat;
  else if (geom_mode == 2 || geom_mode == 3)
    theFoil_mat = composite25pc; // standard 25% each by weight foil
  else
    G4cout << "Invalid geom_mode = " << geom_mode << ", probably about to crash!" << G4endl;

  G4LogicalVolume *theFoil_L = new G4LogicalVolume(theFoil_S, theFoil_mat, "theFoil_L");

  G4VPhysicalVolume *theFoil_P = new G4PVPlacement(0, G4ThreeVector(d_foil, 0, 0), theFoil_L, "theFoil", world_L, false, 0);

  G4VisAttributes *FoilVisAtt = new G4VisAttributes(DARKGRAY);
  FoilVisAtt->SetVisibility(true);
  theFoil_L->SetVisAttributes(FoilVisAtt);


  //----Target slab for validation
  //
  if (geom_mode <= 1) {
    G4double slab_x = 1.0*cm; // D in glass equation
    G4double slab_y = 30.0*cm;
    G4double slab_z = 30.0*cm;

    G4Box *targetSlab_S = new G4Box("targetSlab_S", slab_x/2.0, slab_y/2.0, slab_z/2.0);

    G4Material *targetSlab_mat = (geom_mode == 0 ? U238_simple : NatAl_mat);

    G4LogicalVolume *targetSlab_L = new G4LogicalVolume(targetSlab_S, targetSlab_mat, "targetSlab_L");

    G4VPhysicalVolume *targetSlab_P = new G4PVPlacement(0, G4ThreeVector(50.0*cm, 0, 0), targetSlab_L, "targetSlab", world_L, false, 0);

    targetSlab_L->SetVisAttributes((geom_mode == 0 ? UVisAtt : AlVisAtt));
  }


  //----open-source Black Sea-like object: see, e.g., http://www.pnas.org/content/113/31/8618
  //
  if (geom_mode == 2) {
    G4double r1_Pu_BS = 6.27*cm;
    G4double r2_Pu_BS = 6.70*cm;
    G4double r1_HE_BS = r2_Pu_BS;
    G4double r2_HE_BS = r1_HE_BS + 6.50*cm;
    G4double r1_WU_BS = r2_HE_BS;
    G4double r2_WU_BS = r1_WU_BS + 0.25*cm;
    //G4double r1_Fe_BS = r2_WU_BS;
    //G4double r2_Fe_BS = r1_Fe_BS + 7.00*cm; // ditch the iron

    G4Sphere *Pu_BS_S = new G4Sphere("Pu_BS_S", r1_Pu_BS, r2_Pu_BS, 0, twopi, 0, pi);
    G4Sphere *HE_BS_S = new G4Sphere("HE_BS_S", r1_HE_BS, r2_HE_BS, 0, twopi, 0, pi);
    G4Sphere *WU_BS_S = new G4Sphere("WU_BS_S", r1_WU_BS, r2_WU_BS, 0, twopi, 0, pi);
    //G4Sphere *Fe_BS_S = new G4Sphere("Fe_BS_S", r1_Fe_BS, r2_Fe_BS, 0, twopi, 0, pi);

    G4LogicalVolume *Pu_BS_L = new G4LogicalVolume(Pu_BS_S, WgPu_simple, "Pu_BS_L");
    G4LogicalVolume *HE_BS_L = new G4LogicalVolume(HE_BS_S, HMX,         "HE_BS_L");
    G4LogicalVolume *WU_BS_L = new G4LogicalVolume(WU_BS_S, Uenriched95, "WU_BS_L");
    //G4LogicalVolume *Fe_BS_L = new G4LogicalVolume(Fe_BS_S, NatFe_mat,   "Fe_BS_L");

    G4VPhysicalVolume *Pu_BS_P = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Pu_BS_L, "BlackSeaPlutonium", world_L, false, 0);
    G4VPhysicalVolume *HE_BS_P = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), HE_BS_L, "BlackSeaHighExplosive", world_L, false, 0);
    G4VPhysicalVolume *WU_BS_P = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), WU_BS_L, "BlackSeaUranium", world_L, false, 0);
    //G4VPhysicalVolume *Fe_BS_P = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Fe_BS_L, "BlackSeaIron", world_L, false, 0);

    Pu_BS_L->SetVisAttributes(PuVisAtt);
    HE_BS_L->SetVisAttributes(HEVisAtt);
    WU_BS_L->SetVisAttributes( UVisAtt);
    //Fe_BS_L->SetVisAttributes(AlVisAtt);

    // attach the SDs to their LogicalVolumes
    Pu_BS_L->SetSensitiveDetector(Pu_BS_SD);
    Pu_BS_SD->SetMassOfLogical(Pu_BS_L->GetMass());
    SDMan->AddNewDetector(Pu_BS_SD);

    HE_BS_L->SetSensitiveDetector(HE_BS_SD);
    HE_BS_SD->SetMassOfLogical(HE_BS_L->GetMass());
    SDMan->AddNewDetector(HE_BS_SD);

    WU_BS_L->SetSensitiveDetector(WU_BS_SD);
    WU_BS_SD->SetMassOfLogical(WU_BS_L->GetMass());
    SDMan->AddNewDetector(WU_BS_SD);
  }


  //----mother volume for Black Sea hoax slabs
  //
  if (geom_mode == 3) {
    G4Box *slabV_S = new G4Box("slabV_S", 0.50*m, 0.25*m, 0.25*m);
    G4LogicalVolume *slabV_L = new G4LogicalVolume(slabV_S, Vacuum, "slabV_L");

    G4VisAttributes *SlabVisAtt = new G4VisAttributes();
    SlabVisAtt->SetVisibility(false);
    slabV_L->SetVisAttributes(SlabVisAtt);

    // apply a rotation of theta degrees about the z-axis
    // user can change from 0 to something else
    G4RotationMatrix *slabRotation = new G4RotationMatrix();
    slabRotation->rotateZ(0.0*pi/180.0);

    G4VPhysicalVolume *slabV_P = new G4PVPlacement(slabRotation, G4ThreeVector(0, 0, 0), slabV_L, "slabV", world_L, false, 0);

    //----plates of metal to match the genuine Black Sea object signal -- a geometric hoax detectable by rotation
    //
    double plateY = 35.5*cm;
    double plateZ = 35.5*cm;

    // thicknesses of the plates
    double plate1X = 0.86*cm/2.0; // WgPu_simple
    double plate2X = 13.0*cm/2.0; // HMX
    double plate3X = 0.50*cm/2.0; // Uenriched95
    //double plate4X = 14.0*cm/2.0; // NatFe_mat
    double plate5X = plate1X; // WgPu_simple
    double plate6X = plate2X; // HMX
    double plate7X = plate3X; // Uenriched95
    //double plate8X = plate4X; // NatFe_mat

    // x-locations of the plates
    double plate1d = plate1X/2.0;
    double plate2d = plate1d + plate1X/2.0 + plate2X/2.0;
    double plate3d = plate2d + plate2X/2.0 + plate3X/2.0;
    //double plate4d = plate3d + plate3X/2.0 + plate4X/2.0;

    G4Box *plate1_S = new G4Box("plate1_S", plate1X/2.0, plateY/2.0, plateZ/2.0);
    G4Box *plate2_S = new G4Box("plate2_S", plate2X/2.0, plateY/2.0, plateZ/2.0);
    G4Box *plate3_S = new G4Box("plate3_S", plate3X/2.0, plateY/2.0, plateZ/2.0);
    //G4Box *plate4_S = new G4Box("plate4_S", plate4X/2.0, plateY/2.0, plateZ/2.0);
    G4Box *plate5_S = new G4Box("plate5_S", plate5X/2.0, plateY/2.0, plateZ/2.0);
    G4Box *plate6_S = new G4Box("plate6_S", plate6X/2.0, plateY/2.0, plateZ/2.0);
    G4Box *plate7_S = new G4Box("plate7_S", plate7X/2.0, plateY/2.0, plateZ/2.0);
    //G4Box *plate8_S = new G4Box("plate8_S", plate8X/2.0, plateY/2.0, plateZ/2.0);

    G4LogicalVolume *plate1_L = new G4LogicalVolume(plate1_S, WgPu_simple, "plate1_L");
    G4LogicalVolume *plate2_L = new G4LogicalVolume(plate2_S, HMX, "plate2_L");
    G4LogicalVolume *plate3_L = new G4LogicalVolume(plate3_S, Uenriched95, "plate3_L");
    //G4LogicalVolume *plate4_L = new G4LogicalVolume(plate4_S, NatFe_mat, "plate4_L");
    G4LogicalVolume *plate5_L = new G4LogicalVolume(plate5_S, WgPu_simple, "plate5_L");
    G4LogicalVolume *plate6_L = new G4LogicalVolume(plate6_S, HMX, "plate6_L");
    G4LogicalVolume *plate7_L = new G4LogicalVolume(plate7_S, Uenriched95, "plate7_L");
    //G4LogicalVolume *plate8_L = new G4LogicalVolume(plate8_S, NatFe_mat, "plate8_L");

    G4VPhysicalVolume *plate1_P = new G4PVPlacement(0, G4ThreeVector(plate1d, 0, 0), plate1_L, "plate1", slabV_L, false, 0);
    G4VPhysicalVolume *plate2_P = new G4PVPlacement(0, G4ThreeVector(plate2d, 0, 0), plate2_L, "plate2", slabV_L, false, 0);
    G4VPhysicalVolume *plate3_P = new G4PVPlacement(0, G4ThreeVector(plate3d, 0, 0), plate3_L, "plate3", slabV_L, false, 0);
    //G4VPhysicalVolume *plate4_P = new G4PVPlacement(0, G4ThreeVector(plate4d, 0, 0), plate4_L, "plate4", slabV_L, false, 0);
    G4VPhysicalVolume *plate5_P = new G4PVPlacement(0, G4ThreeVector(-plate1d, 0, 0), plate5_L, "plate5", slabV_L, false, 0);
    G4VPhysicalVolume *plate6_P = new G4PVPlacement(0, G4ThreeVector(-plate2d, 0, 0), plate6_L, "plate6", slabV_L, false, 0);
    G4VPhysicalVolume *plate7_P = new G4PVPlacement(0, G4ThreeVector(-plate3d, 0, 0), plate7_L, "plate7", slabV_L, false, 0);
    //G4VPhysicalVolume *plate8_P = new G4PVPlacement(0, G4ThreeVector(-plate4d, 0, 0), plate8_L, "plate8", slabV_L, false, 0);

    plate1_L->SetVisAttributes(PuVisAtt);
    plate2_L->SetVisAttributes(HEVisAtt);
    plate3_L->SetVisAttributes( UVisAtt);
    //plate4_L->SetVisAttributes(AlVisAtt);
    plate5_L->SetVisAttributes(PuVisAtt);
    plate6_L->SetVisAttributes(HEVisAtt);
    plate7_L->SetVisAttributes( UVisAtt);
    //plate8_L->SetVisAttributes(AlVisAtt);
  }


  //-------------------------------------------------------------------------
  // List of physical volumes

  G4cout << G4endl;
  G4cout << "List of physical volumes: " << G4endl;
  for (int i = 0; i < PVStore->size(); i++) {
    G4VPhysicalVolume *pv = (*PVStore)[i];
    G4cout << "  " << pv->GetName() << " (" << pv->GetLogicalVolume()->GetMaterial()->GetName() << ")" << G4endl;
  }
  G4cout << G4endl;


  return world_P;
}





