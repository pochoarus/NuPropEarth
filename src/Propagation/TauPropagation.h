#ifndef _TAU_PROPAGATION_H_
#define _TAU_PROPAGATION_H_

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHEPEVTParticle.h"

#include "LeptonPropagation.h"

using namespace Tauolapp;

using namespace genie;

namespace genie {

  class TauPropagation : LeptonPropagation {

    public :
      TauPropagation(string proposaltable, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist={});
     ~TauPropagation() {}

      std::vector<GHepParticle> Propagate(GHepParticle * lepton, double minenergy);
      void Propagate(GHepParticle * lepton, double length, double minenergy);

    private :

      std::vector<GHepParticle> Decay(GHepParticle * tau);

      double mass;
      double polarization;

  };

}      // genie namespace
#endif // _TAU_PROPAGATION_H_