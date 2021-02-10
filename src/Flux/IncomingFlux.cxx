#include <cassert>

#include <TMath.h> 
#include <TSystem.h>


#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/PrintUtils.h"

#include "IncomingFlux.h"

using namespace std;
using namespace genie;
using namespace genie::flux;
using namespace genie::constants;

//____________________________________________________________________________
IncomingFlux::IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double cthmono, double emin, double emax, double emono, double radius)
{

  fNewNeutrino = true;

  if (pdg==-1) {
    fPdg.push_back( 12);
    fPdg.push_back(-12);
    fPdg.push_back( 14);
    fPdg.push_back(-14);
    fPdg.push_back( 16);
    fPdg.push_back(-16);
  }
  else fPdg.push_back(pdg);

  fSpectralIndex = alpha;          

  fCThmin = cthmin;
  fCThmax = cthmax;
  fCThmono = cthmono;

  fEmin = emin;          
  fEmax = emax;          
  fEmono = emono;          

  fRadius = radius;

  fNeutrino = new GHepParticle();

  LOG("IncomingFlux", pDEBUG) << "Flux  " << pdg;

}
//___________________________________________________________________________
IncomingFlux::~IncomingFlux()
{

}
//___________________________________________________________________________
PDGCodeList & IncomingFlux::FluxParticles(void)
{

  PDGCodeList * list = new PDGCodeList(false);
  list->push_back( 12);
  list->push_back(-12);
  list->push_back( 14);
  list->push_back(-14);
  list->push_back( 16);
  list->push_back(-16);
  return *list;

}
//___________________________________________________________________________
bool IncomingFlux::GenerateNext(void)
{

  if (fNewNeutrino) this->GenerateNext_1try();
  else fNewNeutrino=true;

  return true;

}
//___________________________________________________________________________
bool IncomingFlux::GenerateNext_1try(void)
{

  RandomGen * rnd = RandomGen::Instance();

  double e = 0.;
  if (fEmono>0.) e = fEmono;
  else {
    if (fSpectralIndex==1) e = exp(log(fEmin)+rnd->RndFlux().Rndm()*log(fEmax/fEmin));
    else {
      double emin = TMath::Power(fEmin,  1.-fSpectralIndex);
      double emax = TMath::Power(fEmax,1.-fSpectralIndex);
      e = TMath::Power(emin+(emax-emin)*rnd->RndFlux().Rndm(),1./(1.-fSpectralIndex));
    }
  }

  double dz = 0.;
  if (fCThmono>0.) dz = fCThmono;
  else             dz = fCThmin+(fCThmax-fCThmin)*rnd->RndFlux().Rndm(); 

  fNeutrino->SetPdgCode(fPdg[rnd->RndFlux().Integer(fPdg.size())]);
  fNeutrino->SetMomentum(0.,e*TMath::Sqrt(1-dz*dz),e*dz,e); 

  double r   = fRadius * TMath::Sqrt( rnd->RndFlux().Rndm() );
  double phi = 2 * kPi * rnd->RndFlux().Rndm();
  double x = r*TMath::Sin(phi);
  double y = r*TMath::Cos(phi);
  fNeutrino->SetPosition( x, y, -fREarth_m, 0. );

  return true;

}
//___________________________________________________________________________
void IncomingFlux::InitNeutrino(GHepParticle Nu)
{

  fNeutrino->Copy(Nu); 
  fNewNeutrino = false; 
      
}