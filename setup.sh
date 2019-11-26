#----------------------------------------------------------
# setup root
#----------------------------------------------------------

export ROOTSYS=/usr/local/root/6.10.02/
source $ROOTSYS/bin/thisroot.sh

#----------------------------------------------------------
# setup everything needed to run/compile PropEarth
# (it's a lot)
#----------------------------------------------------------

#genie                                                                                                                                                                        
export GENIE=/sps/km3net/users/agarcia/software/GENIE-HEDIS-dev/GENIE-HEDIS
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH

#log4cpp (used by genie)                                                                                                               
export LOG4CPPLIB=$KM3NET_THRONG_DIR/src/gSeaGen/log4cpp/lib
export LD_LIBRARY_PATH=$LOG4CPPLIB:$LD_LIBRARY_PATH

#pythia
export PYTHIA6=/pbs/software/centos-7-x86_64/pythia/6.4.28
export LD_LIBRARY_PATH=$PYTHIA6:$LD_LIBRARY_PATH

#lhapdf                                                                                                               
export LHAPDF=/usr/local/lhapdf/6.1.6
export LHAPATH=${GENIE}/data/evgen/pdfs
export PATH=$LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$LHAPDF/lib:$LD_LIBRARY_PATH

#apfel
export APFEL=/sps/km3net/users/agarcia/software/apfel
export LD_LIBRARY_PATH=$APFEL/install/lib:$LD_LIBRARY_PATH

#TAUOLA
export TAUOLA=/sps/km3net/users/agarcia/software/TAUOLA
export LD_LIBRARY_PATH=$TAUOLA/lib:$LD_LIBRARY_PATH

#cern
export CERN=/afs/in2p3.fr/cernlib/amd64_sl7/2006b_oldAFS/x86_64-slc5-gcc43-opt

#nuearthprop
export NUPROPEARTH=/sps/km3net/users/agarcia/software/NuPropEarth
export PATH=$NUPROPEARTH/bin:$PATH
export LD_LIBRARY_PATH=$NUPROPEARTH/lib:$LD_LIBRARY_PATH
