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
IncomingFlux::IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double cthmono, double emin, double emax, double emono)
{

  fNewNeutrino = true;

  fPdg = pdg;          

  fSpectralIndex = alpha;          

  fCThmin = cthmin;
  fCThmax = cthmax;
  fCThmono = cthmono;

  fEmin = emin;          
  fEmax = emax;          
  fEmono = emono;          


  LOG("IncomingFlux", pDEBUG) << "Flux  " << fPdg;

}
//___________________________________________________________________________
IncomingFlux::~IncomingFlux()
{

}
//___________________________________________________________________________
PDGCodeList & IncomingFlux::FluxParticles(void)
{

  PDGCodeList * fList = new PDGCodeList(false);
  fList->push_back(12);
  fList->push_back(-12);
  fList->push_back(14);
  fList->push_back(-14);
  fList->push_back(16);
  fList->push_back(-16);
  return *fList;

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

  // Reset previously generated neutrino code / 4-p / 4-x / wights
  this->ResetSelection();

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

  fgP4I.SetPxPyPzE ( 0., e*TMath::Sqrt(1-dz*dz), e*dz, e );

  return true;

}
//___________________________________________________________________________
void IncomingFlux::InitNeutrino(double px, double py, double pz, double e, double x, double y, double z, double t, int pdg, const TLorentzVector & IntVer)
{

  fgPdgCI = pdg;

  fgX4I.SetXYZT(IntVer.X()+x,IntVer.Y()+y,IntVer.Z()+z,IntVer.T()+t);
  fgP4I.SetPxPyPzE(px,py,pz,e);
  
  fNewNeutrino = false; 
      
}
//___________________________________________________________________________
void IncomingFlux::ResetSelection(void)
{

  fgPdgCI = fPdg;
  fgX4I.SetXYZT    ( 0.,     0., -fREarth_m, 0. );
  
}
