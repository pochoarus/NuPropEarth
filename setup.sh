#----------------------------------------------------------
# setup root
#----------------------------------------------------------

export ROOTSYS=/path/to/your/favourite/ROOT6
source $ROOTSYS/bin/thisroot.sh

#----------------------------------------------------------
# setup genie3 (w/ hedis)
#----------------------------------------------------------

#genie                                                                                                                                                                        
export GENIE=/path/to/directory/were/you/want/to/install/GENIE
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH

#pythia
export PYTHIA6=/path/to/directory/were/PYTHIA6/is/installed
export LD_LIBRARY_PATH=$PYTHIA6:$LD_LIBRARY_PATH

#lhapdf                                                                                                               
export LHAPDF=/path/to/directory/were/LHAPDF6/is/installed
export LHAPATH=${GENIE}/data/evgen/pdfs
export PATH=$LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$LHAPDF/lib:$LD_LIBRARY_PATH

#apfel
export APFEL=/path/to/directory/were/apfel/is/installed
export LD_LIBRARY_PATH=$APFEL/lib:$LD_LIBRARY_PATH

#----------------------------------------------------------
# setup nupropearth
#----------------------------------------------------------

#TAUOLA
export TAUOLA=/path/to/directory/were/tauola/is/installed
export LD_LIBRARY_PATH=$TAUOLA/lib:$LD_LIBRARY_PATH

#cern
export CERN=/path/to/directory/were/cernlib/is/installed

#nuearthprop
export NUPROPEARTH=${PWD}
export PATH=$NUPROPEARTH/bin:$PATH
export LD_LIBRARY_PATH=$NUPROPEARTH/lib:$LD_LIBRARY_PATH
