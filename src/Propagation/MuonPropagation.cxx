
#include "MuonPropagation.h"

MuonPropagation::MuonPropagation(string proposaltable, double ecut, double vcut, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist) : 
    LeptonPropagation(13, proposaltable, ecut, vcut, seed, gd, skiplist) 
{

  SetMass(constants::kMuonMass);

}

std::vector<GHepParticle> MuonPropagation::Propagate(GHepParticle * lepton, double minenergy) 
{ 

  Step(lepton,ComputeDepth(lepton),minenergy);

  std::vector<GHepParticle> PropProd;
  if (lepton->Status()!=kIStDecayedState) PropProd.push_back(*lepton);

  return PropProd;

}

void MuonPropagation::Propagate(GHepParticle * lepton, double length, double minenergy) 
{ 

  Step(lepton,length,minenergy);

  if (lepton->Status()==kIStDecayedState) {
    //decay are not interesting
    lepton->SetPdgCode(0);
  }

}




