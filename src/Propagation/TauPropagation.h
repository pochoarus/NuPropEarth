#ifndef _TAU_PROPAGATION_H_
#define _TAU_PROPAGATION_H_

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/NaturalIsotopes.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include "PROPOSAL/PROPOSAL.h"

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHEPEVTParticle.h"

#include <TSystem.h>

using namespace Tauolapp;

using namespace genie;

struct Ionisation_Constants {
  double I;              ///< ionization potential [eV]
  double C;              ///< ionization formula constants
  double a;
  double m;
  double X0;
  double X1;
  double d0;
};

namespace genie {

  class TauPropagation {

    public :
      TauPropagation(string ptype, int seed, GeomAnalyzerI * gd);
     ~TauPropagation();

      vector<GHepParticle> Propagate(GHepParticle * tau );


    private :

      vector<GHepParticle> Decay(double pdg, double vx, double vy, double vz, double t, double px, double py, double pz, double e);

      void ComputeDepth(GHepParticle * p, double &avgrho, double &lengthi);
      vector<PROPOSAL::Components::Component> GetComponent(map<int,double> composition);
      void ConfigProposal(GeomAnalyzerI * gd);

      RandomGen * rnd;

      GeomAnalyzerI * geom_driver;

      string tauproptype;

      double mtau;
      double d_lifetime;
      double polarization;

      //PROPOSAL
      PROPOSAL::Propagator * ProposalTau;

      //TAUSIC variables
      int TAUIN;
      int TAUMODEL;
      int ITFLAG;

  };

}      // genie namespace
#endif // _TAU_PROPAGATION_H_