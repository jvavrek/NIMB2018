name: README_G4NRF.md  
date: August 14, 2018  
auth: Jayson Vavrek  
mail: jvavrek@mit.edu  
copr: This version (C) 2018 Jayson Vavrek under the MIT License.

G4NRF is a C++ codebase for the simulation of nuclear resonance fluorescence (NRF) physics in Geant4. The package was first developed by David Jordan and Glen Warren at Pacific Northwest National Laboratory, but was largely overhauled by Jayson Vavrek in 2015-2018 (see https://arxiv.org/abs/1807.01701).

For questions/maintenance, submit an issue, submit a pull request, or contact the author directly (email above), though under the MIT License, the author provides no warranty and assumes no liability. This and other copyright notices are subject to change if the codebase is included in the official Geant4 source distribution.

### G4NRF CLASSES ###

* G4NRF: main NRF process class that handles NRF cross section calculation (whether on-the-fly or table-interpolated) and final state generation

* G4NRFPhysics: wrapper class for the NRF process

* G4NRFNuclearLevel: class that holds info about a single nuclear excited state

* G4NRFNuclearLevelManager: class that contains the G4NRFNuclearLevels for a particular isotope

* G4NRFNuclearLevelStore: class that contains all the G4NRFNuclearLevelManagers


### INSTALLATION ###

G4NRF is currently configured to take advantage of the standard Geant4 make process; all the user needs to do is:

1. copy the G4NRF {src,include}/\* files into their respective directories in the user's existing Geant4 application code. This includes not only the class files beginning with G4NRF\*, but also the AngularCorrelation\*, c2_function\*, and RootLinkDef.hh files. See, e.g., the ZKExp directory structure.

2. add the NRF process to the user's existing physicsList:

        #include "G4NRFPhysics.hh"
        ...
        void physicsList::ConstructPhysics() {
        ...
          RegisterPhysics(new G4NRFPhysics("NRF", use_xsec_tables, use_xsec_integration, force_isotropic));
        ...
        }

    (again, see the ZKExp example code). Here the user can either fix the boolean variables use\_xsec\_tables, use\_xsec\_integration, and force\_isotropic, or (as in ZKExp), allow them to be input as command line arguments.

3. unpack the NRF_Database tarball somewhere (if not already done so when following README.md for the ZKExp example code) and add

        export G4NRFGAMMADATA=/absolute/path/to/NRF_Database

    to ~/.bashrc or ~/.bash_profile
