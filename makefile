#############################################################################
#
# name: makefile
# date: August 14, 2018
# auth: Jayson Vavrek, Zach Hartwig
# mail: jvavrek@mit.edu, hartwig@psfc.mit.edu
#
# desc: This file is the GNU makefile that controls the ZKExp
#       build system. Users beware: this is no ordinary Geant4 code!
#       The build system handles a number fancy maneuvers including
#       generating dictionaries for data readout into ROOT framework.
#
# dpnd: 0. The ROOT toolkit     (mandatory)
#       1. Geant4 build with Qt (optional)
#
#############################################################################


##################
#  G4 Makefile  #
##################

# Set the'name' and G4TARGET variables. Note that the 'name' variable is required
# to set the name of the build directory in $G4INSTALL/tmp/$G4SYSTEM directory.
name := ZKExp
ZKLIBS := libZKExpRoot.so
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif


.PHONY: all
all: lib/$(ZKLIBS) lib bin

debug: CXXFLAGS += -Wall -Wextra -pedantic -fno-omit-frame-pointer -O0 -g -ggdb
debug: all

include $(G4INSTALL)/config/binmake.gmk

# Once the G4 Makefile is included, we can customize compilers (if desired)
# CC=clang
# CXX=clang++
CXX=g++

# Newest version of G4 with ROOT results in massive warnings output
# since ROOT local variable 's' shadows the G4Unit 's' for
# seconds. This doesn't affect anything so suppress warning.
CXXFLAGS := $(subst -Wshadow,,$(CXXFLAGS))

# Deal with Geant4 version parsing. Inelegant, but should work
G4V := $(shell geant4-config --version)
G4VS := $(subst ., ,$(G4V))
G4MAJV := $(word 1,$(G4VS))
G4MINV := $(word 2,$(G4VS))
G4SUBV := $(word 3,$(G4VS))
CXXFLAGS += -DG4MAJV=$(G4MAJV) -DG4MINV=$(G4MINV) -DG4SUBV=$(G4SUBV)


##########
#  ROOT  #
##########

# ROOT classes are presently used in ZKExp for their immense
# utility; this requires compiling and linking against ROOT as a
# dependency. Use 'root-config' to obtain the header and library dirs
ROOTINCLUDES = -I$(shell root-config --incdir)
ROOTDISTLIBS = $(shell root-config --nonew --libs --glibs)

CPPFLAGS += $(ROOTINCLUDES)
LDLIBS += $(ROOTDISTLIBS) -L./lib -lZKExpRoot

ZK_ROOT_FILES = $(wildcard include/*.rhh)

lib/libZKExpRoot.so : lib/ZKExpDict.o 
	@echo -e "\nBuilding the ZKExp ROOT library ...\n"
	$(CXX) $(ROOTDISTLIBS) -shared -o $@  $^

lib/ZKExpDict.o : lib/ZKExpDict.cc
	@echo -e "Building $@ ..."
	$(CXX) $(CXXFLAGS) $(ROOTINCLUDES) -I. -c -o $@ $<

lib/ZKExpDict.cc : $(ZK_ROOT_FILES) include/RootLinkDef.hh
	@echo -e "Generating the ZKExp ROOT dictionary ..."
	rootcint -f $@ -c $^


.PHONY:

libclean:
	@echo -e "\nCleaning up the ZKExp libraries ...\n"
	@rm -f lib/*
