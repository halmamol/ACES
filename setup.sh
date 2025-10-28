#!/bin/bash

cd /mnt/netapp2/Store_uni/home/usc/ie/dcr/software/ROOT/ROOT_6.26.14/install/bin
source thisroot.sh

cd /mnt/netapp2/Store_uni/home/usc/ie/dcr/software/GEANT4/geant4.10.3.3/install/bin
source geant4.sh

cd /mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/WCSim/build/mydir
source bin/this_wcsim.sh

cd /mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/ACES
export BONSAIDIR=/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/hk-BONSAI
export ROOTSYS=/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/ROOT/ROOT_6.26.14/install
export WCSIM_BUILD_DIR=/mnt/netapp2/Store_uni/home/usc/ie/dcr/software/hk/WCSim/build/mydir