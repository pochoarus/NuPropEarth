#ifndef _MUON_PROPAGATION_H_
#define _MUON_PROPAGATION_H_

#include "LeptonPropagation.h"

using namespace genie;

namespace genie {

  class MuonPropagation : LeptonPropagation {

    public :
      MuonPropagation(string proposaltable, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist={});
     ~MuonPropagation() {}

      std::vector<GHepParticle> Propagate(GHepParticle * lepton, double minenergy);
      void Propagate(GHepParticle * lepton, double length, double minenergy);

  };

}      // genie namespace
#endif // _MUON_PROPAGATION_H_