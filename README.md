# NuPropEarth

NuPropEarth is a tool used to propagated neutrinos through the Earth. It has the structure of a general-purpose Monte Carlo event generator, and therefore allows following the path and interactions of individual neutrinos as they travel through Earth on an an event-by-event basis. The main goal of the NuPropEarth framework is to compute the attenuation coefficients, defined as ratios between the incoming neutrino flux to the Earth and the flux arriving at the detector volume.

## Authors: 

<pre>
Alfonso Garcia < alfonsog \at nikhef.nl >
Rhorry Gauld < r.gauld \at nikhef.nl >
Aart Heijboer < aart.heijboer \at nikhef.nl >
Juan Rojo < j.rojo \at vu.nl >
</pre>


## Prerequisites

These are the main dependencies:

- GENIE3 (w/ HEDIS)
- TAUOLA++ (v1.1.8) [modify src/tauolaCInterfaces/TauolaParticle.cxx line 295 using double m=Tauola::getTauMass();]
- TAUSIC

Then GENIE3 (w/ HEDIS) requires several external packages:

- ROOT6
- PYTHIA6
- LHAPDF6
- APFEL


## INSTALL

1. Install all the external packages: PYTHIA6,TAUOLA,TAUSIC,LHAPDF6,APFEL,ROOT6

2. Download NuPropEarth module on your machine 

```
git clone https://github.com/pochoarus/NuPropEarth.git
cd NuPropEarth
```

3. Define the enviroment in which you will work in source.sh (for instance GENIE,PYHTIA6,LHAPDF6,etc.)

4. Source the enviroment

```
source setup.sh
```

5. Install this version of [GENIE3 (w/ HEDIS)](https://github.com/pochoarus/GENIE-HEDIS/tree/nupropearth)

```
git clone -b nupropearth https://github.com/pochoarus/GENIE-HEDIS.git $GENIE
cd $GENIE
./configure --enable-lhapdf6 --enable-apfel --with-lhapdf6-inc=$LHAPDF/include --with-lhapdf6-lib=$LHAPDF/lib --with-apfel-inc=$APFEL/include --with-apfel-lib=$APFEL/lib
make
```

6. Install NuPropEarth in your machine

```
cd $NUPROPEARTH
make
```


## EXAMPLE

To compute the attenuation you just have to run

```
NUMBEROFEVENTS="1e3" #to have enough statistic use 1e6
NUPDG="14" #14=numu, 12=nue, 16=nutau, -14=anumu, -12=anue, -16=anutau, 
COSTHETA="0.1"
TUNE="GHE19_00a_00_000" #GHE19_00a_00_000=BGR18(member=0), GHE19_00b_00_000=CSMS11(member=0)
cd $NUPROPEARTH
source setup.sh
ComputeAttenuation -n $NUMBEROFEVENTS -t $COSTHETA -p $NUPDG --cross-sections ${GENIE}/genie_xsec/${TUNE}_dx0.01dy0.01_n200_nuall_nucleon_noglres.xml --event-generator-list HEDIS --tune $TUNE
```


## Citation

To cite this work, and for more information, please refer to

Complete predictions for high-energy neutrino propagation in matter

preprint: [arXiv:2004.04756](https://arxiv.org/abs/2004.04756)

NuPropEarth is GENIE based application [GENIE citing rules](https://hep.ph.liv.ac.uk/~costasa/genie/citing.html).



