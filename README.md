![logo](/logo.png)


NuPropEarth is a tool used to propagated neutrinos through the Earth. It has the structure of a general-purpose Monte Carlo event generator, and therefore allows following the path and interactions of individual neutrinos as they travel through Earth on an an event-by-event basis. The main goal of the NuPropEarth framework is to compute the attenuation coefficients, defined as ratios between the incoming neutrino flux to the Earth and the flux arriving at the detector volume.

## Authors: 

<pre>
Alfonso Garcia < alfonsog \at nikhef.nl >
Rhorry Gauld < r.gauld \at nikhef.nl >
Aart Heijboer < aart.heijboer \at nikhef.nl >
Juan Rojo < j.rojo \at vu.nl >
Victor Valera < vvalera \at nbi.ku.dk >
</pre>


## Prerequisites

These are the main dependencies:

- GENIE (3.2.0 or above)
- TAUOLA++ (v1.1.8)
- PROPOSAL (v6)

Then GENIE requires several external packages:

- ROOT6
- PYTHIA6
- LHAPDF6


## Installation instructions

Tested in CentOS Linux release 7.9.2009 (Core) with access to sft.cern.ch (CVFMS).

```
source /cvmfs/sft.cern.ch/lcg/releases/gcc/7.3.0-cb1ee/x86_64-centos7/setup.sh

export WORKDIR=/path/to/installation/directory/
mkdir $WORKDIR

#install tauola
cd $WORKDIR
wget https://tauolapp.web.cern.ch/tauolapp/resources/TAUOLA.1.1.8/TAUOLA.1.1.8.tar.gz
tar xf TAUOLA.1.1.8.tar.gz
rm TAUOLA.1.1.8.tar.gz
mv TAUOLA tauola
cd tauola
sed -i "s/double m=.*/double m=Tauola::getTauMass();/g" src/tauolaCInterfaces/TauolaParticle.cxx
./configure --without-hepmc --without-hepmc3
make

#install proposal
cd $WORKDIR
git clone --recursive -b 6.1.6 https://github.com/tudo-astroparticlephysics/PROPOSAL proposal
cd proposal/
mkdir build install
cd build
export PATH=/cvmfs/sft.cern.ch/lcg/releases/CMake/3.8.2-ece19/x86_64-centos7-gcc7-opt/bin:$PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/CMake/3.8.2-ece19/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-std=c++11 -DADD_PYTHON=OFF  -DCMAKE_INSTALL_PREFIX=../install
cmake --build . --target install

#install genie
cd $WORKDIR
git clone -b R-3_04_00 https://github.com/GENIE-MC/Generator.git genie
cd genie/
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-abd9a/x86_64-centos7-gcc7-opt/bin/thisroot.sh
export GENIE=$PWD
./configure --disable-lhapdf5 --enable-lhapdf6 --with-lhapdf6-inc=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-4b9c6/x86_64-centos7-gcc7-opt/include --with-lhapdf6-lib=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-4b9c6/x86_64-centos7-gcc7-opt/lib --with-pythia6-lib=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia6/429.2-63d8b/x86_64-centos7-gcc7-opt/lib --with-libxml2-inc=/cvmfs/sft.cern.ch/lcg/releases/libxml2/2.9.7-830a9/x86_64-centos7-gcc7-opt/include/libxml2 --with-libxml2-lib=/cvmfs/sft.cern.ch/lcg/releases/libxml2/2.9.7-830a9/x86_64-centos7-gcc7-opt/lib --with-log4cpp-inc=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/log4cpp/2.8.3-aeffd/x86_64-centos7-gcc7-opt/include --with-log4cpp-lib=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/log4cpp/2.8.3-aeffd/x86_64-centos7-gcc7-opt/lib
sed -i "s/-lPythia6/-lpythia6/g" src/make/Make.include
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia6/429.2-63d8b/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/tbb/2018_U1-d3621/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/GSL/2.1-36ee5/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/log4cpp/2.8.3-aeffd/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/libxml2/2.9.7-830a9/x86_64-centos7-gcc7-opt/lib:$LD_LIBRARY_PATH
make -j16


#install nupropearth
cd $WORKDIR
git clone https://github.com/pochoarus/NuPropEarth.git nupropearth
cd nupropearth/
export NUPROPEARTH=$PWD
export PROPOSAL=${WORKDIR}/proposal/install
export TAUOLA=${WORKDIR}/tauola
make -j16


# Creating enviroment script
cat > ${WORKDIR}/setup.sh << EOL
source /cvmfs/sft.cern.ch/lcg/releases/gcc/7.3.0-cb1ee/x86_64-centos7/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.12.04-abd9a/x86_64-centos7-gcc7-opt/bin/thisroot.sh
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia6/429.2-63d8b/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-4b9c6/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/log4cpp/2.8.3-aeffd/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/tbb/2018_U1-d3621/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/GSL/2.1-36ee5/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/libxml2/2.9.7-830a9/x86_64-centos7-gcc7-opt/lib:\$LD_LIBRARY_PATH
export GENIE=$GENIE
export LD_LIBRARY_PATH=\$GENIE/lib:\$LD_LIBRARY_PATH
export PATH=\$GENIE/bin:\$PATH
export LHAPATH=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.2.1-4b9c6/x86_64-centos7-gcc7-opt/share/LHAPDF/:/cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/pdfs/:\$GENIE/data/evgen/pdfs/
export HEDIS_SF_DATA_PATH=/cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/hedis-sf
export PHOTON_SF_DATA_PATH=/cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/photon-sf
export XSEC_SPLINES=/cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/splines
export TAUOLA=$TAUOLA
export LD_LIBRARY_PATH=\$TAUOLA/lib:$LD_LIBRARY_PATH
export PROPOSAL=$PROPOSAL
export LD_LIBRARY_PATH=\$PROPOSAL/lib64:$LD_LIBRARY_PATH
export NUPROPEARTH=$NUPROPEARTH
export LD_LIBRARY_PATH=\$NUPROPEARTH/lib:\$LD_LIBRARY_PATH
export PATH=\$NUPROPEARTH/bin:\$PATH
export GXMLPATH=\$NUPROPEARTH
EOL
```



## Example

Generate muon neutrino interactions in oxygen using a power law spectrum
```
source $WORKDIR/setup.sh

#first we create geometry file
cat > ./oxygen_sphere.xml << EOL
<?xml version="1.0" encoding="ISO-8859-1"?>
<earth_prem>
  <param_set name="LayerA" radius="100" density="1.,0.,0.,0.">
    <param element="1000080160">   1.0           </param>
  </param_set>
</earth_prem>
EOL

BuildEarth -xml ./oxygen_sphere.xml -root ./oxygen_sphere.root

#generate events
TUNE=GHE19_00a_00_000 #bgr model
VertexGenerator --seed 1 --output ./test.root --number-of-events 1e3 --probe 14 --alpha 1 --costheta -1 --energy 1e2,1e10 --offset 0 --detector-radius 1000 --detector-height 1000 --detector-position "0,0,0" --geometry-limit 100e3 --geometry ./oxygen_sphere.root --event-generator-list CCHEDIS --tune $TUNE --cross-sections ${XSEC_SPLINES}/${TUNE}.xml
```

Propagate neutrinos through Earth (PREM model)
```
```

## Citation

To cite this work, and for more information, please refer to

Complete predictions for high-energy neutrino propagation in matter

preprint: [arXiv:2004.04756](https://arxiv.org/abs/2004.04756)

NuPropEarth is GENIE based application [GENIE citing rules](https://hep.ph.liv.ac.uk/~costasa/genie/citing.html).



