name: README.md  
date: August 14, 2018  
auth: Jayson Vavrek, Zach Hartwig  
mail: jvavrek at mit dot edu, hartwig at psfc dot mit dot edu  
copr: This version (C) 2018 Jayson Vavrek under the MIT License.

The purpose of ZKExp is to provide the Geant4 simulation infrastructure
for the MIT LNSP Physical Cryptography (formerly "Zero Knowledge") project:

1. http://www.pnas.org/content/113/31/8618
2. https://arxiv.org/abs/1712.02904
3. https://arxiv.org/abs/1807.01701

Moreover, this release uses ZKExp as a publicly-available demonstration of
the updated G4NRF code base, which simulates nuclear resonance fluorescence
(NRF) in Geant4. See also README_G4NRF.md for further detail.

ZKExp is based on a slick Geant4 user application template provide by Zach
Hartwig. In addition to typical Geant4 geometry, physics (including G4NRF),
particle source, and user action classes, ZKExp includes a particularly
key feature: the ability to incorporate ROOT classes seamlessly (including
automatic ROOT dictionary generation) into the user's Geant4 code. This provides,
for example, the ability to readout data into TTrees, save all data persistently
for later analysis in TFiles, create histograms on the fly, and so on.


### DIRECTORY STRUCTURE

* ZKExp.cc : standard G4 top-level file containing mandatory C++ main() function.

* makefile : GNU makefile for building ZKExp

* README.md : this file

* README_G4NRF.md : documentation for the G4NRF package in particular

* include/ : standard class header files (*.hh) including G4NRF headers

* src/ : standard class source files (*.cc) including G4NRF classes

* lib/ : contains the ROOT library and dictionary files

* runtime/ : contains run and visualization macros

* NRF_Database.tar.gz : gzip-ed tarball of NRF level and gamma data (from ENSDF)

* brems_distributions.root : ROOT input file for the bremsstrahlung and sampling
                             spectra for the PGA

* setup.sh : setup script for environmental variables

* createSampling.C : root script that creates brem_normalized.root and an
                     importance-sampled distribution in hSample.root;
                     requires converter\_2p7MeV\_50M\_bremTotal\_option4.root
                     in order to be re-run, but this file is too large to
                     include in the GitHub repository.


### INSTALLATION ###

The user should first configure their environment and compile
ZKExp to ensure everything is working correctly. To do so:

1. Ensure your system contains an up-to-date, complete installation of
   Geant4 (http://geant4.cern.ch/) with the extended libraries. It is
   optional to build the Qt libraries and include them when building
   the Geant4 libraries. Ensure that your environment and G4 working
   directory paths are correctly configured.

2. Ensure your system contains an up-to-date, complete installation of
   the ROOT data analysis toolkit (http://root.cern.ch/drupal/) and that
   your environment has been correctly configured. NOTE: ZKExp has only
   been minimally tested with ROOT 6!

3. Clone or download ZKExp from its github repository

4. Unpack the NRF_Database tarball to a suitable location

5. Add the following lines to your .bashrc or .bash_profile file:

        export ZKEXP_TOPDIR=/absolute/path/to/ZKExp
        source $ZKEXP_TOPDIR/setup.sh
        export G4NRFGAMMADATA=/absolute/path/to/NRF_Database

6. Restart the terminal or re-source .bashrc / .bash_profile


### BUILDING ###

Building the ZKExp simulation is presently handled by a GNU
makefile located in the top-level directory.

To build ZKExp simply type:
       
    $ make -jX

where 'X' is number of cores to use in the compilation process.

Typing 'make clean' will clean up all the transient build files while
'make libclean' will clean up the ROOT library and dictionary build
files contained in the lib/ directory. Typing 'make debug' will build
the executable with several debug flags allowing the user to use a
debugger such as GDB or a profiler such as Mac Instruments.


### EXECUTION ###

To run the ZKExp executable, type:

    $ ZKExp <vis> <macro> <seed> <geom> <beam> <xsec_tables> <xsec_integration> <force_iso>

where:

* vis   = {visOn, visQt, visOff} to select visualization option
* macro = path to runtime macro to run immediately  
* seed  = integer: random seed AND unique ROOT output file identifier  
* geom  = integer: {0: U-238 obj and foil, 1: same but Al-27, 2: Black Sea obj and composite foil, 3: same but slab geom}  
* beam  = integer: {0: U-238 resonance sampling, 1: Al-27 resonance sampling, 2: input histo res. samp.}  
* xsec\_tables = {0, 1}: true/false whether to use NRF cross section tables instead of on-the-fly evaluation  
* xsec\_integration = {0, 1}: true/false whether to use numerically-integrated NRF cross section (else Gaus approx)  
* force\_iso = {0, 1}: true/false whether to force isotropic angular correlation for NRF emission  


### EXAMPLES ###

Some example ZKExp configurations, including those used in Paper \#3 above:

    $ ZKExp visOn

starts the Geant4 UI with OpenGL visualization and no run macro.

    $ ZKExp visOff runtime/ZK.mac 31415 0 0 1 1 0

runs the U-238 geometry with appropriate beam using tabulated numerically-integrated NRF cross sections with isotropic emission.

    $ ZKExp visOff runtime/ZK.mac 31415 1 1 1 1 0

will run the analogous Al-27 simulation.

    $ ZKExp visOff runtime/ZK.mac 31415 3 2 1 1 1

will run the 'Black Sea hoax' simulation with a simulated bremsstrahlung beam and angular correlations activated.

NOTE: there are no flags at this time to change target thicknesses. The user will need to edit the source code (src/geometryConstruction.cc) directly.


### OUTPUT ###

ZKExp uses the ROOT TTree class to record information about photons stepping out of the G4LogicalVolume 'theFoil\_L', which is the reference/encryption foil in a transmission NRF experiment. Cuts of theta > 3pi/4 and E > 1.7 MeV are enforced (but can be changed in src/steppingAction.cc). Information about surviving photons is written to the TFile ZKExp\_outfile\_\<seed\>.root, or to ZK\_\<timestamp\>.root if no seed is specified. The user can write energy information from ROOTTreeFoil to a TH1D energy spectrum via, e.g.,

    $ root ZKExp_outfile_31415.root
    root> TH1D *h = new TH1D("h", "h", 1000, 1.7, 2.7);
    root> h->Sumw2();
    root> ROOTTreeFoil->Draw("energy>>h","weightF","e")

The 'weightF' variable is used to keep track of statistical weights arising from importance sampling (only with beam = 2). In this case, it is important for the user to explicity create a TH1D histogram instead of letting TTree::Draw() automatically create a TH1F, since double precision will prevent roundoff error that can easily be introduced by floating point arithmetic when the weights span several orders of magnitude.

     
### ADAPTATION ###

In order to use ZKExp as a starting point for a user's
simulation, he/she will need to take a few minor steps. Let's assume
the user wants to develop a simulation named
"ZKPrototype". He/she should take the following steps:

1. Copy ZKExp directory to a directory named "ZKPrototype" and cd
   into the newly named directory.

2. Change "ZKExp.cc" to "ZKPrototype.cc"

3. Replace all occurrences of "ZKExp" with "ZKPrototype" in the
   following files:
  * makefile
  * ZKPrototype.cc

4. Replace all occurrences of "ZKEXP" with "ZKPROTOTYPE" in
   "setup.sh" file

5. Update the .bashrc file with the environmental variables

6. Restart terminal / resource .bashrc for changes to take effect
