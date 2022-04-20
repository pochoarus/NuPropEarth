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
IncomingFlux::IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double emin, double emax, double detpos[3], double radius, double height, double offset, bool towards)
{

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

  fEmin = emin;          
  fEmax = emax;          
  fDetPos[0]  = detpos[0];
  fDetPos[1]  = detpos[1];
  fDetPos[2]  = detpos[2];
  fDetRadius  = radius;
  fDetHeight  = height;         
  fOffset     = offset;         
  fTowardsDet = towards;         

  fNeutrino = new GHepParticle();

  fGlobalWeight = fPdg.size();
  fGlobalWeight *= (fCThmin==fCThmax) ? 1. : 2.*TMath::Pi()*(fCThmax-fCThmin);
  if (fEmin==fEmax) fGlobalWeight *= 1.;
  else              fGlobalWeight *= (fSpectralIndex==1) ? TMath::Log(fEmax/fEmin) : (TMath::Power(fEmax,(1.-fSpectralIndex))-TMath::Power(fEmin,(1.-fSpectralIndex)))/(1.-fSpectralIndex);

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
void IncomingFlux::InitNeutrino(void)
{

  RandomGen * rnd = RandomGen::Instance();

  double e = 0.;
  if (fEmin==fEmax) e = fEmin;
  else {
    if (fSpectralIndex==1) e = TMath::Exp(TMath::Log(fEmin)+rnd->RndFlux().Rndm()*TMath::Log(fEmax/fEmin));
    else {
      double emin = TMath::Power(fEmin,  1.-fSpectralIndex);
      double emax = TMath::Power(fEmax,  1.-fSpectralIndex);
      e = TMath::Power(emin+(emax-emin)*rnd->RndFlux().Rndm(),1./(1.-fSpectralIndex));
    }
  }

  fNeutrino->SetPdgCode(fPdg[rnd->RndFlux().Integer(fPdg.size())]);

  double dz    = (fCThmin==fCThmax) ? fCThmin : fCThmin+(fCThmax-fCThmin)*rnd->RndFlux().Rndm();
  double theta = TMath::ACos(dz);
  double phi   = 2. * TMath::Pi() * rnd->RndFlux().Rndm();
  double dx    = TMath::Sin(theta) * TMath::Sin(phi);
  double dy    = TMath::Sin(theta) * TMath::Cos(phi);

  double ShiftPos[3] = { 0., 0., 0. };

  if ( fDetRadius!=0 || fDetHeight!=0 ) {
    
    // If the detector has a finite extension we use a cylindrical shape, depending on the angle we choose a surface to generate the initial position (top, bottom, or lateral face)
    double A1 = TMath::Pi()*TMath::Power(fDetRadius,2.)*TMath::Abs(dz); // Corresponds to the area of a circular shape projected into a plane at an angle given by fCth
    double A2 = 2.0*fDetRadius*fDetHeight*TMath::Sqrt(1.-dz*dz); // Corresponds to the area of the lateral face of the cylinder projected into the plane

    double f = A1/(A1+A2);  // This fraction is used to select an starting surface by comparing it with a random number between 0 and 1
  
    if (rnd->RndFlux().Rndm() < f)  { // Event at the top/bottom of the cylinder
      double randR = fDetRadius*TMath::Sqrt(rnd->RndFlux().Rndm());
      double randAlpha = 2.*TMath::Pi()*rnd->RndFlux().Rndm();
      ShiftPos[2] = (dz < 0) ? -fDetHeight*0.5 : fDetHeight*0.5; 
      ShiftPos[1] = randR*TMath::Sin(randAlpha);
      ShiftPos[0] = randR*TMath::Cos(randAlpha);
    }
    else {
      ShiftPos[2] = fDetHeight*(rnd->RndFlux().Rndm()-0.5);
      ShiftPos[1] = 2.*fDetRadius*(rnd->RndFlux().Rndm()-0.5);
      ShiftPos[0] = -TMath::Sqrt(fDetRadius*fDetRadius - ShiftPos[1]*ShiftPos[1]);
      TVector3 vaux(ShiftPos[0],ShiftPos[1],ShiftPos[2]);
      vaux.RotateZ( TMath::ATan2( -dy, -dx ) );
      ShiftPos[2] = vaux.Z();
      ShiftPos[1] = vaux.Y();
      ShiftPos[0] = vaux.X();
    }

    LOG("IncomingFlux", pDEBUG) << "Shift: " << ShiftPos[0] << "  " << ShiftPos[1] << "  " << ShiftPos[2];

    fAgen = A1+A2;

  } 

  fNeutrino->SetPosition( fDetPos[0]+ShiftPos[0]+fOffset*dx,  fDetPos[1]+ShiftPos[1]+fOffset*dy, fDetPos[2]+ShiftPos[2]+fOffset*dz, 0. );

  if (fTowardsDet) fNeutrino->SetMomentum( -e*dx, -e*dy, -e*dz, e );
  else             fNeutrino->SetMomentum(  e*dx,  e*dy,  e*dz, e );

}
//___________________________________________________________________________
void IncomingFlux::InitNeutrino(GHepParticle Nu)
{

  fNeutrino->Copy(Nu); 
        
}
