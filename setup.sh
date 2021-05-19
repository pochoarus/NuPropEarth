#----------------------------------------------------------
# setup root
#----------------------------------------------------------

export ROOTSYS=/path/to/directory/where/ROOT6/is/installed
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

#----------------------------------------------------------
# setup genie3 (w/ hedis)
#----------------------------------------------------------

#genie                                                                                                                                                                        
export GENIE=/path/to/directory/where/GENIE/is/installed
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE/lib:$LD_LIBRARY_PATH

#log4cpp
export LOG4CPP=/path/to/directory/where/LOG4CPP/is/installed
export LD_LIBRARY_PATH=$LOG4CPP/lib:$LD_LIBRARY_PATH

#libxml2
export LIBXML2=/path/to/directory/where/LIBXML2/is/installed
export LD_LIBRARY_PATH=$LIBXML2/lib:$LD_LIBRARY_PATH

#pythia
export PYTHIA6=/path/to/directory/where/PYTHIA6/is/installed
export LD_LIBRARY_PATH=$PYTHIA6:$LD_LIBRARY_PATH

#lhapdf                                                                                                               
export LHAPDF=/path/to/directory/where/LHAPDF6/is/installed
export LHAPATH=${GENIE}/data/evgen/pdfs
export PATH=$LHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$LHAPDF/lib:$LD_LIBRARY_PATH

#apfel (not mandatory)
export APFEL=/path/to/directory/where/apfel/is/installed
export LD_LIBRARY_PATH=$APFEL/lib:$LD_LIBRARY_PATH

#----------------------------------------------------------
# setup nupropearth
#----------------------------------------------------------

#TAUOLA
export TAUOLA=/path/to/directory/where/tauola/is/installed
export LD_LIBRARY_PATH=$TAUOLA/lib:$LD_LIBRARY_PATH

#TAUSIC
export TAUSIC=/path/to/directory/where/tausic/is/installed

#cern
export CERN_LEVEL=/version/of/cernlib
export CERN=/path/to/directory/where/cernlib/is/installed
export CERN_ROOT=$CERN/$CERN_LEVEL
export PATH=$CERN_ROOT/bin:$PATH

#nuearthprop
export NUPROPEARTH=/path/to/directory/where/nupropearth/is/installed
export PATH=$NUPROPEARTH/bin:$PATH
export LD_LIBRARY_PATH=$NUPROPEARTH/lib:$LD_LIBRARY_PATH
export GXMLPATH=$NUPROPEARTH

