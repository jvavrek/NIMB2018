/////////////////////////////////////////////////////////////////////////////////
//
// name: ZKExp.cc
// date: August 14, 2018
// auth: Jayson Vavrek, Zach Hartwig
// mail: jvavrek@mit.edu, hartwig@psfc.mit.edu
// copr: This version (C) 2018 Jayson Vavrek under the MIT License.
//
// desc: The purpose of ZKExp is to provide a basic but mature
//       Geant4 code structure and build system to function as the
//       basis for creating all Geant4 simulations associated with the
//       MIT Zero Knowledge (ZK) experiment. In addition to typical
//       Geant4 geometry, physics, particle source, and user action
//       classes, ZKExp includes two particularly key features:
//
//       0. Ability to incorporate ROOT classes seamlessly (including
//          automatic ROOT dictionary generation).
//
//       1. Inclusion of a singleton manager class to handle all ROOT
//          input/ouput of important simulation data.
//
//      A working ROOT 5 installation on the system is mandatory.
//      ZKExp has not been tested with ROOT 6!
//
//      In addition, cmd line options are available to specify Qt,
//      OpenGL, or no visualization options, and runtime files placed
//      in the runtime/ directory are automatically available from
//      with a running instance of the ZKExp simulation
//
//      It is intended that the user copy this entire directory to a
//      new directory, replace "ZKExp" with his/her own
//      simulation name, and use this as a foundation to build his/her
//      simulation.
//
/////////////////////////////////////////////////////////////////////////////////

// C classes
#include <unistd.h>
#include <sys/resource.h>
#include <stdlib.h>

// C++ classes
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>

// Geant4 classes
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIcsh.hh"
#include "G4SDManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// ZK Classes
#include "geometryConstruction.hh"
#include "physicsList.hh"
#include "PGA.hh"
#include "steppingAction.hh"
#include "stackingAction.hh"
#include "eventAction.hh"
#include "runAction.hh"
#include "rootStorageManager.hh"


int main(int argc, char *argv[]) {
  if (getenv("ZKEXP_TOPDIR") == NULL) {
    G4cout << "\nZK ANNOUNCEMENT: The environmental variable ZKEXP_TOPDIR must be set to\n"
     <<   "                 the top-level ZKExp directory in the current shell session.\n"
     <<   "                 It is recommended to source the ZKExp setup file 'setup.sh' in\n"
     <<   "                 your .bashrc file!\n"
     << G4endl;
    G4Exception("ZK main()", "ZKMainException001", FatalException, "ZK EXCEPTION: Evasive maneuvers!\n");
  }

  G4String ZK_TOPDIR = (G4String)std::getenv("ZKEXP_TOPDIR");
  G4String ZK_RUNTIMEDIR = ZK_TOPDIR + "/runtime";

  ////////////////////////////////
  // Parse command line options //
  ////////////////////////////////

  // Set default settings before going into the argument parsing
  G4bool visualization = false;
  G4bool visQt = false;

  G4String macro_file = "";

  G4long new_seed = 0;
  G4int geom_mode = 0;
  G4int beam_mode = 0;

  G4bool use_xsec_tables = true;
  G4bool use_xsec_integration = true;
  G4bool force_isotropic = false;

  // Template command-line call:
  // $ ZKExp <vis> <macro> <seed> <geom> <beam> <xsec_tables> <xsec_integration> <force_iso>
  //
  // where:
  //  - vis   = {visOn, visQt, visOff} to select visualization
  //  - macro = path to runtime macro to run immediately
  //  - seed  = integer: random seed AND unique ROOT output file identifier
  //  - geom  = integer: {0: U-238 obj and foil, 1: same but Al-27, 2: Black Sea obj and composite foil, 3: same but slab geom}
  //  - beam  = integer: {0: U-238 resonance sampling, 1: Al-27 resonance sampling, 2: input histo res. samp.}
  //  - xsec_tables = {0, 1}: true/false whether to use NRF cross section tables instead of on-the-fly evaluation
  //  - xsec_integration = {0, 1}: true/false whether to use numerically-integrated NRF cross section (else Gaus approx)
  //  - force_iso = {0, 1}: true/false whether to force isotropic angular correlation for NRF emission
  //
  // Example:
  // $ ZKExp visOff runtime/ZK.mac 31415 0 0 1 1 0
  // runs the DU geometry with appropriate beam using tabulated but numerically-integrated NRF cross sections

  if (argc > 1) {
    // visualization
    if (strcmp(argv[1], "visOn") == 0) {
      visualization = true;
    } else if (strcmp(argv[1], "visQt") == 0) {
      visualization = true;
      visQt = true;
    }
  }
  if (argc > 2) { macro_file = argv[2]; }                         // macro file
  if (argc > 3) { new_seed  = (G4long) atol(argv[3]); }           // random seed
  if (argc > 4) { geom_mode = atoi(argv[4]); }                    // geometry mode
  if (argc > 5) { beam_mode = atoi(argv[5]); }                    // beam mode
  if (argc > 6) { use_xsec_tables = (atoi(argv[6]) == 1); }       // build cross section tables?
  if (argc > 7) { use_xsec_integration = (atoi(argv[7]) == 1); }  // use numerically-integrated NRF xsec?
  if (argc > 8) { force_isotropic = (atoi(argv[8]) == 1); }       // force isotropic NRF angular correlation?

  // bounds checking on new_seed, beam_mode is in PGA, while geom_mode is in geometryConstruction

  //////////////////////////////////////////////
  // Initialize mandatory/user Geant4 classes //
  //////////////////////////////////////////////

  // Create the runManager to handle program flow
  G4RunManager *runManager = new G4RunManager;

  // Assign the mandatory user-derived classes to the runManager and
  // initialize it before creation of the user actions so that it can
  // be subsequently accessed from the constructors of the user action
  // classes.
  geometryConstruction *theGC = new geometryConstruction(geom_mode);
  runManager->SetUserInitialization(theGC);

  physicsList *thePL = new physicsList(use_xsec_tables, use_xsec_integration, force_isotropic);
  runManager->SetUserInitialization(thePL);

  PGA *thePGA = new PGA(new_seed, beam_mode);
  runManager->SetUserAction(thePGA);

  runManager->Initialize();

  // Create the user action classess and assign to the run manager
  steppingAction *SteppingAction = new steppingAction();
  SteppingAction->SetBeamDir(thePGA->GetBeamDir());
  runManager->SetUserAction(SteppingAction);

  stackingAction *StackingAction = new stackingAction();
  StackingAction->SetBeamDir(thePGA->GetBeamDir());
  runManager->SetUserAction(StackingAction);

  eventAction *EventAction = new eventAction();
  runManager->SetUserAction(EventAction);

  runAction *RunAction = new runAction();
  runManager->SetUserAction(RunAction);


  /////////////////////////////////////////////////////////////////
  // Initialize the user-interface manager with defaults/aliases //
  /////////////////////////////////////////////////////////////////

  // Get the pointer to the U(ser) I(nterface) manager.
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  {
    G4String alias = "/control/alias";
    G4String exe = "/control/execute {run}/";

    // Set the alias to contain absolute path to runtime files
    UImanager->ApplyCommand(alias + " run " + ZK_RUNTIMEDIR);
    UImanager->ApplyCommand(alias + " mac " + exe + "ZK.mac");

    // Set some basic Geant4-specific verbosity defaults
    UImanager->ApplyCommand("/run/verbose 1");
    UImanager->ApplyCommand("/event/verbose 0");
    UImanager->ApplyCommand("/hits/verbose 0");
    UImanager->ApplyCommand("/tracking/verbose 0");
    UImanager->ApplyCommand("/control/verbose 1");
  }


  /////////////////////////////////
  // Initialize the ROOT manager //
  /////////////////////////////////

  rootStorageManager *theRSManager = new rootStorageManager(new_seed);


  // ZKExp presently always randomizes its RNG seed; changing
  // conditional to false will start all simulations with same seed
  if (true)
    CLHEP::HepRandom::setTheSeed(time(0) + 42);

  // Only include visualization capabilities if Geant4 has been
  // built with the G4VIS_USE flag
#ifdef G4VIS_USE
  G4VisManager *visManager = new G4VisExecutive((visualization ? "warnings" : "quiet"));
  visManager->Initialize();

  // Color particle trajectories by particle type
  G4TrajectoryDrawByParticleID *colorModel = new G4TrajectoryDrawByParticleID;
  colorModel->Set("neutron", "cyan");
  colorModel->Set("gamma", "green");
  colorModel->Set("e-", "red");
  colorModel->Set("e+", "blue");
  colorModel->Set("proton", "yellow");
  colorModel->SetDefault("gray");
  visManager->RegisterModel(colorModel);
  visManager->SelectTrajectoryModel(colorModel->Name());
#endif

  // At present, two options are provided for visualization and
  // control of ZK. First, the hotness of G4UI graphical user
  // interface using Qt; second standard command-line-plus-OpenGL.
  if (visQt) {
    // Build the G4UI GUI and run the Qt visualization macro
#ifdef G4UI_USE
    G4UIExecutive *UIexecutive = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute {run}/ZK.Qt.vis");
    UImanager->ApplyCommand("/vis/scene/add/trajectories");
    UImanager->ApplyCommand("/vis/scene/add/hits");
    UIexecutive->SessionStart();
    delete UIexecutive;
#endif
  } else {
    if (visualization) {
      // Run the OpenGL visualization macro
      UImanager->ApplyCommand("/control/execute {run}/ZK.OGLIX.vis");
      UImanager->ApplyCommand("/vis/scene/add/trajectories");
      UImanager->ApplyCommand("/vis/scene/add/hits");
    }

    // Create a decent 'tcsh'-like prompt for tab completion, command
    // history, etc.  Also, style points for cooler prompt
    G4String prompt = "ZK >> ";
    G4int maxHist = 200;
    G4UIsession* session = new G4UIterminal(new G4UItcsh(prompt, maxHist));


    // As Gallagher said: "Styyyyyyyyyyyle!"
    G4cout << "\n\n"
           << "\t\t    ZZZZZZZ    K    K             \n"
           << "\t\t          Z    K   K              \n"
           << "\t\t         Z     K  K               \n"
           << "\t\t        Z      KKK                \n"
           << "\t\t     ZZZZZ     K  K               \n"
           << "\t\t      Z        K   K              \n"
           << "\t\t     Z  ERO    K    K NOWLEDGE    \n"
           << "\t\t    Z          K    K             \n"
           << "\t\t    ZZZZZZ     K    K             \n"
           << "\n\n      *******    WELCOME TO THE ZKExp SIMULATION    *******\n\n\n"
           << G4endl;

    // If a macro is specified (e.g., in batch mode), skip the terminal prompt and apply
    // the macro immediately
    if (macro_file == "") {
      session->SessionStart();
    } else {
      G4String macroCmd = "/control/execute " + macro_file;
      UImanager->ApplyCommand(macroCmd);
    }
    delete session;

#ifdef G4VIS_USE
  if (visualization)
    delete visManager;
#endif
  }

  // General garbage collection
  delete runManager;
  delete theRSManager;
  G4cout << "\n" << G4endl;

  return 0;
}

