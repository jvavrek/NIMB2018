#!/bin/bash
#
# name: setup.sh
# date: 17 Sep 14
# auth: Zach Hartwig
# 
# desc: Bash script to configure the user's environment for running
#       the ZK Geant4 simulation
#
export PATH=$ZKEXP_TOPDIR/scripts:$PATH
export LD_LIBRARY_PATH=$ZKEXP_TOPDIR/lib:$LD_LIBRARY_PATH
