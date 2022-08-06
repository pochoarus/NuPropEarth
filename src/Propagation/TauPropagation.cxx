
#include "TauPropagation.h"

RandomGen * taurnd;

double TauRandomGenerator() { return taurnd->RndGen().Rndm(); }

TauPropagation::TauPropagation(string proposaltable, double ecut, double vcut, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist) :
    LeptonPropagation(15, proposaltable, ecut, vcut, seed, gd, skiplist)
{

  LOG("TauPropagation", pDEBUG) << "Initializing Random...";
  taurnd = RandomGen::Instance();

  LOG("TauPropagation", pDEBUG) << "Initializing TAUOLA...";
  Tauola::setSeed(seed,0,0);
  Tauola::setRandomGenerator(TauRandomGenerator);
  Tauola::initialize();

  SetMass(Tauola::getTauMass());
  polarization = 1; //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped

}

std::vector<GHepParticle> TauPropagation::Propagate(GHepParticle * lepton, double minenergy) { 

  double length = ComputeDepth(lepton);

  Step(lepton,length,minenergy);

  std::vector<GHepParticle> PropProd;
  if (lepton->Status()==kIStDecayedState) PropProd = Decay(lepton);
  else                                    PropProd.push_back(*lepton);

  return PropProd;

}

void TauPropagation::Propagate(GHepParticle * lepton, double length, double minenergy) { 

  Step(lepton,length,minenergy);

  if (lepton->Status()==kIStDecayedState) {
    //look for muon
    std::vector<GHepParticle> PropProd = Decay( lepton );
    for ( auto & tpr : PropProd ) {
        if ( pdg::IsMuon(TMath::Abs(tpr.Pdg())) ) { lepton->Copy(tpr); break; }
        //if muon not found then we dont care about this product
        lepton->SetPdgCode(0);
    }
  }

}


std::vector<GHepParticle> TauPropagation::Decay(GHepParticle * tau) {

  //decay tau
  TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
  TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( tau->Pdg(), 1, tau->Px(), tau->Py(), tau->Pz(), tau->E(), Tauola::getTauMass(), -1, -1, -1, -1 );
  Tauola_evt->addParticle(Tauola_tau);
  Tauola::decayOne(Tauola_tau,true,0.,0.,polarization);

  std::vector<GHepParticle> DecProd;
  
  for ( int sec=1; sec<Tauola_evt->getParticleCount(); sec++ ) {
    int spdg   = Tauola_evt->getParticle(sec)->getPdgID();
    double se  = Tauola_evt->getParticle(sec)->getE();
    double spx = Tauola_evt->getParticle(sec)->getPx();
    double spy = Tauola_evt->getParticle(sec)->getPy();
    double spz = Tauola_evt->getParticle(sec)->getPz();
    LOG("TauPropagation", pDEBUG) << "Product: " << spdg << ", E = " << se << " GeV";
    LOG("TauPropagation", pDEBUG) << "  Position   = [ " << tau->Vx() << " m, " << tau->Vy() << " m, " << tau->Vz() << " m, " << tau->Vt() << " s ]";
    LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << spx/se << " , " << spy/se << " , " << spz/se << " ]";
    DecProd.push_back(GHepParticle(spdg,kIStUndefined,-1,-1,-1,-1,spx,spy,spz,se,tau->Vx(),tau->Vy(),tau->Vz(),tau->Vt()));
  }

  delete Tauola_evt;

  return DecProd;

}





