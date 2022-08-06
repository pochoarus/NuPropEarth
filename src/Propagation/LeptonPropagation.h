#ifndef _LEPTON_PROPAGATION_H_
#define _LEPTON_PROPAGATION_H_

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/NaturalIsotopes.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include "PROPOSAL/PROPOSAL.h"

#include <TSystem.h>

using namespace std;
using namespace genie;
using namespace genie::geometry;

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

  class LeptonPropagation {

    public :
      virtual ~LeptonPropagation() {}

      virtual std::vector<GHepParticle> Propagate(GHepParticle * lepton, double minenergy) = 0;
      virtual void Propagate(GHepParticle * lepton, double length, double minenergy) = 0;

      std::vector<GHepParticle> StepShowers(GHepParticle * lepton, double length, double minenergy); 

    protected :

      LeptonPropagation(int pdg, string proposaltable, double ecut, double vcut, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist);

      double ComputeDepth(GHepParticle * p);

      void SetMass(double m) { mass = m; } 

      void Step(GHepParticle * lepton, double length, double minenergy); 
      

    private :

      std::vector<PROPOSAL::Components::Component> GetComponent(map<int,double> composition);
      void ConfigProposal(int pdg, string proposaltable, double ecut, double vcut, vector<string> skiplist);

      ROOTGeomAnalyzer * geom_driver;

      PROPOSAL::Propagator * Proposal;

      double mass;

  };

}      // genie namespace
#endif // _LEPTON_PROPAGATION_H_