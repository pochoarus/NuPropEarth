#ifndef _TAU_PROPAGATION_H_
#define _TAU_PROPAGATION_H_

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHEPEVTParticle.h"

#include <TSystem.h>

using namespace Tauolapp;

using namespace genie;

namespace genie {

  class TauPropagation {

    public :
      TauPropagation(string ptype, int seed, GeomAnalyzerI * gd);
     ~TauPropagation();

      vector<GHepParticle> Propagate(GHepParticle * tau );
      vector<GHepParticle> Decay(double pdg, double vx, double vy, double vz, double t, double px, double py, double pz, double e);

      void ComputeDepth(GHepParticle * p, double &avgrho, double &lengthi);
      void SetRockDensity(double rd) { rockdensity = rd; }

    private :

      RandomGen * rnd;

      GeomAnalyzerI * geom_driver;

      string tauproptype;

      double mtau;
      double d_lifetime;
      double polarization;

      //TAUSIC variables
      double rockdensity;
      int TAUIN;
      int TAUMODEL;
      int ITFLAG;

  };

}      // genie namespace
#endif // _TAU_PROPAGATION_H_