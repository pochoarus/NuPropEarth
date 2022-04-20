#ifndef _HADRON_PROPAGATION_H_
#define _HADRON_PROPAGATION_H_

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include <TPythia6.h>
#include <TMCParticle.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>


using namespace genie;
using namespace genie::geometry;

namespace genie {

  class HadronPropagation {

    public :
      HadronPropagation(ROOTGeomAnalyzer * gd);
     ~HadronPropagation(){}

      std::vector<GHepParticle> Propagate(GHepParticle * hadron, string pclass);


    private :

      std::vector<GHepParticle> Decay(double pdg, double kf, double vx, double vy, double vz, double t, double dx, double dy, double dz, double e);

      double HadronInteraction(string pclass, double mass, double E, double rho);
      double HadronInelasticity(string pclass);

      double ComputeDepth(GHepParticle * p);
      double ComputeRho(GHepParticle * p);

      RandomGen * rnd;

      ROOTGeomAnalyzer * geom_driver;

      TPythia6     * fPythia;
      TDatabasePDG * PdgDB;


  };

}      // genie namespace
#endif // _HADRON_PROPAGATION_H_