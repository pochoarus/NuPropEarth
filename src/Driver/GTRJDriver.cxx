#include <cassert>

#include <TVector3.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/GEVGPool.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"

#include "GTRJDriver.h"

using namespace genie;
using namespace genie::constants;

double c = kLightSpeed/(units::meter/units::second);

//____________________________________________________________________________
GTRJDriver::GTRJDriver()
{

  fEventGenList       = "Default";  // <-- set of event generators to be loaded by this driver

  fUnphysEventMask = new TBits(GHepFlags::NFlags()); //<-- unphysical event mask
  for(unsigned int i = 0; i < GHepFlags::NFlags(); i++) fUnphysEventMask->SetBitNumber(i, true);

  fFluxDriver         = 0;     // <-- flux driver
  fGeomAnalyzer       = 0;     // <-- geometry driver
  fGPool              = 0;     // <-- pool of GEVGDriver event generation drivers

  fCurTgtPdg          = 0;
  fCurEvt             = 0;
  fCurdL              = 0.;

  // Force early initialization of singleton objects that are typically
  // would be initialized at their first use later on.
  // This is purely cosmetic and I do it to send the banner and some prolific
  // initialization printout at the top.
  assert( Messenger::Instance()     );
  assert( AlgConfigPool::Instance() );

}
//___________________________________________________________________________
GTRJDriver::~GTRJDriver()
{

  if (fUnphysEventMask) delete fUnphysEventMask;
  if (fGPool) delete fGPool;

}
//___________________________________________________________________________
void GTRJDriver::SetEventGeneratorList(string listname)
{
  LOG("GTRJDriver", pNOTICE)
       << "Setting event generator list: " << listname;

  fEventGenList = listname;
}
//___________________________________________________________________________
void GTRJDriver::UseFluxDriver(GFluxI * flux_driver)
{
  fFluxDriver = flux_driver;
}
//___________________________________________________________________________
void GTRJDriver::UseGeomAnalyzer(ROOTGeomAnalyzer * geom_analyzer)
{
  fGeomAnalyzer = geom_analyzer;
}
//___________________________________________________________________________
void GTRJDriver::Configure(double emin, double emax)
{
  LOG("GTRJDriver", pNOTICE) << utils::print::PrintFramedMesg("Configuring GTRJDriver");

  if (fGPool) delete fGPool;
  fGPool = new GEVGPool;

  // Get the list of flux neutrinos from the flux driver
  LOG("GTRJDriver", pNOTICE) << "Asking the flux driver for its list of neutrinos";
  PDGCodeList fNuList = fFluxDriver->FluxParticles();
  LOG("GTRJDriver", pNOTICE) << "Flux particles: " << fNuList;

  // Get the list of target materials from the geometry driver
  LOG("GTRJDriver", pNOTICE) << "Asking the geometry driver for its list of targets";
  PDGCodeList fTgtList = fGeomAnalyzer->ListOfTargetNuclei();
  LOG("GTRJDriver", pNOTICE) << "Target materials: " << fTgtList;

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;
  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
    for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {

      int target_pdgc   = *tgtiter;
      int neutrino_pdgc = *nuiter;

      InitialState init_state(target_pdgc, neutrino_pdgc);

      LOG("GTRJDriver", pNOTICE)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

      GEVGDriver * evgdriver = new GEVGDriver;
      evgdriver->SetEventGeneratorList(fEventGenList); // specify list of generators
      evgdriver->Configure(init_state);

      Range1D_t rE = evgdriver->ValidEnergyRange();
      if ( emin<rE.min || emax>rE.max ) {
        LOG("GTRJDriver", pFATAL) << "Valid Energy range: " << rE.min << "  " << rE.max;
        LOG("GTRJDriver", pFATAL) << "Define Energy range: " << emin << "  " << emax;
        exit(1);
      }
      evgdriver->CreateSplines( 1000, emax, true ); 
      evgdriver->CreateXSecSumSpline( 1000, emin, emax, true );

      LOG("GTRJDriver", pDEBUG) << "Adding new GEVGDriver object to GEVGPool";
      fGPool->insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );
    } // targets
  } // neutrinos

  LOG("GTRJDriver", pNOTICE) << "Finished configuring GTRJDriver\n\n";

}
//___________________________________________________________________________
int GTRJDriver::GenerateEvent(void)
{
// attempt generating a neutrino interaction by firing a single flux neutrino
//

  fCurEvt    = 0;
  fCurTgtPdg = 0;
  fCurdL     = 0.;

  // Generate a neutrino using the input GFluxI & get current pdgc/p4/x4
  if(!fFluxDriver->GenerateNext()) {
     LOG("GTRJDriver", pFATAL)
         << "*** The flux driver couldn't generate a flux neutrino!!";
     exit(1);
  }

  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();
  TVector3 udir = nup4.Vect().Unit();

  std::vector< pair<double, const TGeoMaterial*> > MatLengths = fGeomAnalyzer->ComputeMatLengths(nux4,nup4);

  double length = 0;
  for ( auto sitr : MatLengths ) length += sitr.first;
  length -= 10; //safety factor of 10 m

  double xf = nux4.X()+length*udir.X();
  double yf = nux4.Y()+length*udir.Y();
  double zf = nux4.Z()+length*udir.Z();
  if ( fGeomAnalyzer->GetGeometry()->FindNode(xf,yf,zf)->GetNumber()==1 ) {    
    LOG("GTRJDriver", pWARN) << "** Final position could be outside the geometry.";
    LOG("GTRJDriver", pWARN) << "** It can happen with events very close to boundaries.";
    LOG("GTRJDriver", pWARN) << "E = " << nup4.E();
    LOG("GTRJDriver", pWARN) << "Pos   = [ " << nux4.X() << " m, " << nux4.Y() << " m, " << nux4.Z() << " m]";
    LOG("GTRJDriver", pWARN) << "Dir   = [ " << udir.X() << ", " << udir.Y() << ", " << udir.Z() << " ]";
    LOG("GTRJDriver", pWARN) << "Final = [ " << xf << " m, " << yf << " m, " << zf << " m]";
    return 0;
  }

  LOG("GTRJDriver", pNOTICE)
     << "\n [-] Generated flux neutrino: "
     << "\n  |----o PDG-code   : " << fFluxDriver->PdgCode()
     << "\n  |----o 4-momentum : " << utils::print::P4AsString(&fFluxDriver->Momentum())
     << "\n  |----o 4-position : " << utils::print::X4AsString(&fFluxDriver->Position());

  // Check if the neutrino interacts in the geometry
  if ( !this->ComputeInteraction( fFluxDriver->PdgCode(), nup4.E(), MatLengths ) ) {
    LOG("GTRJDriver", pNOTICE) << "** Neutrino didnt interact in the volume";
    return 2;
  }

  // Ask the GEVGDriver object to select and generate an interaction and
  // its kinematics for the selected initial state & neutrino 4-momentum
  this->GenerateEventKinematics();
  if(!fCurEvt) {
    LOG("GTRJDriver", pWARN) << "** Couldn't generate kinematics for selected interaction";
    LOG("GTRJDriver", pWARN) << "E = " << nup4.E();
    LOG("GTRJDriver", pWARN) << "Pos   = [ " << nux4.X() << " m, " << nux4.Y() << " m, " << nux4.Z() << " m]";
    LOG("GTRJDriver", pWARN) << "Dir   = [ " << udir.X() << ", " << udir.Y() << ", " << udir.Z() << " ]";
    return 0;
  }

  return 1;

}
//___________________________________________________________________________
bool GTRJDriver::ComputeInteraction(int nupdg, double Enu, std::vector< pair<double, const TGeoMaterial*> > MatLengths)
{

  RandomGen * rnd = RandomGen::Instance();
  double R = rnd->RndEvg().Rndm();
  LOG("GTRJDriver", pDEBUG) << "Rndm [0,1] = " << R;
  double N = -TMath::Log(R);
  LOG("GTRJDriver", pDEBUG) << "N = " << N;

  fCurdL = 0.;
  double NIntCum=0;
  
  for ( auto sitr : MatLengths ) {

    double length            = sitr.first;
    const  TGeoMaterial* mat = sitr.second;
    double rho               = mat->GetDensity() * (units::m3/units::cm3) ;

    LOG("GTRJDriver", pDEBUG) << mat->GetName() << " -> rho = " << rho << " gr/m3";

    double xsec = 0;
    if (mat->IsMixture()) {
      const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
      for (int i = 0; i < mixt->GetNelements(); i++) {
        int mpdg = fGeomAnalyzer->GetTargetPdgCode(mixt, i);
        double Frac = mixt->GetWmixt()[i]; // relative proportion by mass
        xsec += Frac*this->EvalXsec(mpdg,nupdg,Enu);
        LOG("GTRJDriver", pDEBUG) << "tgt: " << mpdg << ", frac =  " << Frac;
      }
    }
    else {
      int mpdg = fGeomAnalyzer->GetTargetPdgCode(mat);
      LOG("GTRJDriver", pDEBUG) << "tgt: " << mpdg;
      xsec = this->EvalXsec(mpdg,nupdg,Enu);
    }

    LOG("GTRJDriver", pDEBUG) << "length: " << length << " m";
    LOG("GTRJDriver", pDEBUG) << "xsec: " << xsec/units::m2 << " m2";

    double NInt  = rho * kNA * xsec/units::m2 * length; // number of interactions in the path

    LOG("GTRJDriver", pDEBUG) << "NInt: " << NInt;

    NIntCum += NInt;

    LOG("GTRJDriver", pDEBUG) << "NIntCum: " << NIntCum;

    if (N>NIntCum) fCurdL += length; //no interaction in this layer -> add length
    else {

      if (mat->IsMixture()) {
        double RFrac = rnd->RndEvg().Rndm();
        LOG("GTRJDriver", pDEBUG) << "RFrac: " << RFrac;
        double FracCum = 0.;
        const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
        for (int i = 0; i < mixt->GetNelements(); i++) {
          FracCum += mixt->GetWmixt()[i]; // relative proportion by mass
          if (RFrac<FracCum) {
            fCurTgtPdg = fGeomAnalyzer->GetTargetPdgCode(mixt, i);
            break;
          }
        }
      }
      else fCurTgtPdg = fGeomAnalyzer->GetTargetPdgCode(mat);

			if(fCurTgtPdg==0) {
        LOG("GTRJDriver", pERROR) << "** Rejecting current flux neutrino (failed to select tgt!)";
        return false;
  		}

      fCurdL += length * (N - (NIntCum-NInt)) / NInt;

      LOG("GTRJDriver", pNOTICE) << "Interaction happening in material = " << fCurTgtPdg;
      LOG("GTRJDriver", pNOTICE) << "Distance to origin = " << fCurdL << " m";

      return true;

    }

  }

  return false;

}
//___________________________________________________________________________
void GTRJDriver::GenerateEventKinematics(void)
{

  int                    nupdg = fFluxDriver->PdgCode();
  const TLorentzVector & nup4  = fFluxDriver->Momentum();
  const TLorentzVector & nux4  = fFluxDriver->Position ();

  // Find the GEVGDriver object that generates interactions for the
  // given initial state (neutrino + target)
  InitialState init_state(fCurTgtPdg, nupdg);
  GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
  if(!evgdriver) {
     LOG("GTRJDriver", pFATAL)
       << "No GEVGDriver object for init state: " << init_state.AsString();
     exit(1);
  }

  // propagate current unphysical event mask 
  evgdriver->SetUnphysEventMask(*fUnphysEventMask);

  // Ask the GEVGDriver object to select and generate an interaction for
  // the selected initial state & neutrino 4-momentum
  LOG("GTRJDriver", pNOTICE)
          << "Asking the selected GEVGDriver object to generate an event";
  fCurEvt = evgdriver->GenerateEvent(nup4);

  // Generate an 'interaction position' in the selected material (in the
  // detector coord system), along the direction of nup4 & set it 
  TVector3 udir = nup4.Vect().Unit();

  TLorentzVector Vtx(nux4.X()+fCurdL*udir.X(), nux4.Y()+fCurdL*udir.Y(), nux4.Z()+fCurdL*udir.Z(), nux4.T() + fCurdL/c);

  fCurEvt->SetVertex(Vtx);

}
//___________________________________________________________________________
double GTRJDriver::EvalXsec(int mpdg, int nupdg, double Enu)
{

    // find the GEVGDriver object that is handling the current init state
    InitialState init_state(mpdg, nupdg);
    GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
    if(!evgdriver) {
      LOG("GTRJDriver", pFATAL)
       << "\n * The MC Job driver isn't properly configured!"
       << "\n * No event generation driver could be found for init state: " 
       << init_state.AsString();
      exit(1);
    }

    const Spline * totxsecspl = evgdriver->XSecSumSpline();
    if(!totxsecspl) {
        LOG("GTRJDriver", pFATAL)
          << "\n * The MC Job driver isn't properly configured!"
          << "\n * Couldn't retrieve total cross section spline for init state: " 
          << init_state.AsString();
        exit(1);
    }

    return totxsecspl->Evaluate( Enu )/pdg::IonPdgCodeToA(mpdg);

}
