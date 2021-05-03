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
IncomingFlux::IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double cthmono, double emin, double emax, double emono, double depth, double radius, double height)
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
  fDepth = depth;
  fRadius = radius;
  fHeight = height;         

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
  if (fCThmono>-2.) dz = fCThmono;
  else             dz = fCThmin+(fCThmax-fCThmin)*rnd->RndFlux().Rndm(); 


  fNeutrino->SetPdgCode(fPdg[rnd->RndFlux().Integer(fPdg.size())]);

  fNeutrino->SetMomentum( 0., e*TMath::Sqrt(1-dz*dz), e*dz, e );

  if (fRadius == 0 && fHeight ==0) {
    fNeutrino->SetPosition( 0.,  0., -fREarth_m + fDepth, 0. );
  }
  else {
    // If the detector has a finite extension we use a cylindrical shape, depending on the angle we choose a surface to generate the initial position (top, bottom, or lateral face)
    double A1 = TMath::Pi()*TMath::Power(fRadius,2.)*TMath::Abs(dz); // Corresponds to the area of a circular shape projected into a plane at an angle given by fCth
    double A2 = 2.0*fRadius*fHeight*TMath::Sqrt(1-dz*dz); // Corresponds to the area of the lateral face of the cylinder projected into the plane

    double f = A1/(A1+A2);  // This fraction is used to select an starting surface by comparing it with a random number between 0 and 1
  
    if (rnd->RndFlux().Rndm() < f)  {
      double randR = fRadius*TMath::Sqrt(rnd->RndFlux().Rndm());
      double randAlpha = 2.*TMath::Pi()*rnd->RndFlux().Rndm();
    
      double randX = randR*TMath::Cos(randAlpha);
      double randY = randR*TMath::Sin(randAlpha);
      
      if (dz < 0) {
        fgX4I.SetXYZT    ( randX,  randY, -fREarth_m + fDepth - fHeight*0.5, 0. ); // Event at the top of the cylinder
      }
      else {
        fgX4I.SetXYZT    ( randX,  randY, -fREarth_m + fDepth + fHeight*0.5, 0. ); // Event on the bottom of the cylinder
      }  
    }
    else {
      double randX = 2.*fRadius*(rnd->RndFlux().Rndm()-0.5);
      double randY = TMath::Sqrt(fRadius*fRadius - randX*randX);
      fNeutrino->SetPosition( randX,  randY, -fREarth_m + fDepth + fHeight*(rnd->RndFlux().Rndm()-0.5), 0. ); // Event at the lateral face
    }
  } 
  
  return true;

}
//___________________________________________________________________________
void IncomingFlux::InitNeutrino(GHepParticle Nu)
{

  fNeutrino->Copy(Nu); 
  
  fNewNeutrino = false; 
      
}
