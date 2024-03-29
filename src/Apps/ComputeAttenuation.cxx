#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TDatabasePDG.h>

#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include "Driver/GTRJDriver.h"
#include "Flux/IncomingFlux.h"
#include "Propagation/TauPropagation.h"
#include "Propagation/HadronPropagation.h"


using std::string;
using std::pair;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

//functions
void GetCommandLineArgs (int argc, char ** argv);
void FillParticle       (GHepParticle * part, int &pdg, double &px, double &py,double &pz, double &e, double &vx, double &vy, double &vz) {
  pdg = part->Pdg();
  px  = part->Px(); py = part->Py(); pz = part->Pz(); e = part->E();
  vx  = part->Vx(); vy = part->Vy(); vz = part->Vz();
}

// User-specified options:
string          gOptOutName     = "./test.root";               
string          gOptGeometry    = "";               
double          gOptGeoLimit    = 0.; //Sphere where all interaction must be contained (usually the maximum radius of geometry)
int             gOptRanSeed     = 0;               
string          gOptInpXSecFile = "";           
int             gOptNev         = 0;              
int             gOptPdg         = -1; //all flavors                   
double          gOptEmin        = 1e2;
double          gOptEmax        = 1e10;
double          gOptCthmin      = -1.;
double          gOptCthmax      =  1.;
double          gOptAlpha       = 1; //flat spectrum in log10e
double          gOptDetPos[3]   = { 0., 0., 0. }; //Detector position at center of the volume
double          gOptDetRadius   = 0.; //Detector Radius (radius of a cylindrical shape) originally a point
double          gOptDetHeight   = 0.; //Detector height (height of the cylinder centerd at gOptDetPos[2]) originally a point
double          gOptOffset      = 0.; //Move the initial positon of neutrino this lenght
string          gOptTauProp     = ""; //Path to proposal tables
bool            gOptHadronProp  = true;

const int NMAXTRIALS = 10;
const Range1D_t spline_Erange = { 1e2, 1e12 }; //energy range limit based on HEDIS splines 

//**************************************************************************
//**************************************************************************
//THIS IS THE MAIN FUNCTION
//**************************************************************************
//**************************************************************************
int main(int argc, char** argv)
{

  // Parse command line arguments
  LOG("ComputeAttenuation", pDEBUG) << "Reading options...";
  GetCommandLineArgs(argc,argv);

  geometry::ROOTGeomAnalyzer * rgeom = new geometry::ROOTGeomAnalyzer(gOptGeometry);
  rgeom -> SetLengthUnits  ( genie::utils::units::UnitFromString("m")     );
  rgeom -> SetDensityUnits ( genie::utils::units::UnitFromString("g_cm3") );
  rgeom -> SetTopVolName   ("");

  TGeoVolume * topvol = rgeom->GetGeometry()->GetTopVolume();
   if(!topvol) {
    LOG("ComputeAttenuation", pFATAL) << " ** Null top ROOT geometry volume!";
    exit(1);
  }
  LOG("ComputeAttenuation", pNOTICE) << topvol->GetName();

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Tau Propagation...";
  TauPropagation * tauprop = new TauPropagation(gOptTauProp,-1,0.001,gOptRanSeed,rgeom,{"MatVacuum"}); //skip vacum from propsal configuration

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Hadron Propagation...";
  HadronPropagation * hadronprop = new HadronPropagation(rgeom);
  TDatabasePDG * PdgDB  = TDatabasePDG::Instance();

  LOG("ComputeAttenuation", pDEBUG) << "Creating GFluxI...";
  IncomingFlux * flx_driver = new IncomingFlux(gOptPdg, gOptAlpha, gOptCthmin, gOptCthmax, gOptEmin, gOptEmax, gOptDetPos, gOptDetRadius, gOptDetHeight, gOptOffset, false);

  RunOpt::Instance()->BuildTune();

  // Iinitialization of random number generators, cross-section table, messenger, cache etc...
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  LOG("ComputeAttenuation", pDEBUG) << "Creating GTRJDriver...";
  GTRJDriver * trj_driver = new GTRJDriver;
  trj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  trj_driver->UseFluxDriver(dynamic_cast<GFluxI *>(flx_driver));
  trj_driver->UseGeomAnalyzer(rgeom);

  LOG("ComputeAttenuation", pDEBUG) << "Configuring GTRJDriver...";
  trj_driver->Configure(spline_Erange.min,spline_Erange.max);

  // create file so tree is saved inside
  TFile * outfile = new TFile(gOptOutName.c_str(),"RECREATE");
  
  int    NInt;
  int    In_Pdg;
  double In_X,In_Y,In_Z;
  double In_PX,In_PY,In_PZ,In_E;
  int    Out_Pdg;
  double Out_X,Out_Y,Out_Z;
  double Out_PX,Out_PY,Out_PZ,Out_E;

  TTree *influx = new TTree("InFlux","InFlux");
  influx->Branch("In_Pdg", &In_Pdg, "In_Pdg/I" );       
  influx->Branch("In_X",   &In_X,    "In_X/D"  );       
  influx->Branch("In_Y",   &In_Y,    "In_Y/D"  );       
  influx->Branch("In_Z",   &In_Z,    "In_Z/D"  );       
  influx->Branch("In_PX",  &In_PX,   "In_PX/D" );       
  influx->Branch("In_PY",  &In_PY,   "In_PY/D" );       
  influx->Branch("In_PZ",  &In_PZ,   "In_PZ/D" );       
  influx->Branch("In_E",   &In_E,    "In_E/D"  );       

  TTree *outflux = new TTree("OutFlux","OutFlux");
  outflux->Branch("NInt",    &NInt,    "NInt/I"    );       
  outflux->Branch("In_Pdg",  &In_Pdg,  "In_Pdg/I"  );       
  outflux->Branch("In_X",    &In_X,    "In_X/D"    );       
  outflux->Branch("In_Y",    &In_Y,    "In_Y/D"    );       
  outflux->Branch("In_Z",    &In_Z,    "In_Z/D"    );       
  outflux->Branch("In_PX",   &In_PX,   "In_PX/D"   );       
  outflux->Branch("In_PY",   &In_PY,   "In_PY/D"   );       
  outflux->Branch("In_PZ",   &In_PZ,   "In_PZ/D"   );       
  outflux->Branch("In_E",    &In_E,    "In_E/D"    );       
  outflux->Branch("Out_Pdg", &Out_Pdg, "Out_Pdg/I" );       
  outflux->Branch("Out_X",   &Out_X,   "Out_X/D"   );       
  outflux->Branch("Out_Y",   &Out_Y,   "Out_Y/D"   );       
  outflux->Branch("Out_Z",   &Out_Z,   "Out_Z/D"   );       
  outflux->Branch("Out_PX",  &Out_PX,  "Out_PX/D"  );       
  outflux->Branch("Out_PY",  &Out_PY,  "Out_PY/D"  );       
  outflux->Branch("Out_PZ",  &Out_PZ,  "Out_PZ/D"  );       
  outflux->Branch("Out_E",   &Out_E,   "Out_E/D"   );       

  LOG("ComputeAttenuation", pDEBUG) << "Event loop...";

  std::vector<GHepParticle> SecNu;
  std::vector<GHepParticle> OutPart;

  int trial = 0; //number of tries for same event
  int NNu = 0;
  NInt = 0;
  while ( NNu<gOptNev ) {

    if ( SecNu.size()==0 ) {
      if ( NNu%100==0 ) LOG("ComputeAttenuation", pNOTICE) << "Event " << NNu << " out of " << gOptNev;
      flx_driver->InitNeutrino();
      FillParticle(flx_driver->GetNeutrino(),In_Pdg,In_PX,In_PY,In_PZ,In_E,In_X,In_Y,In_Z);
    }
    else {
      LOG("ComputeAttenuation", pDEBUG) << "@@SecNu = "   << SecNu[0].Pdg() << " , E = " << SecNu[0].E() << " GeV";
      flx_driver->InitNeutrino(SecNu[0]);
      SecNu.erase(SecNu.begin()); //remove first entry because it is just process
    }

    TVector3 fdir = flx_driver->Momentum().Vect().Unit();
    LOG("ComputeAttenuation", pDEBUG) << "Shooting Neutrino: " << NNu << "( " << NInt << " ) ----->" ;
    LOG("ComputeAttenuation", pDEBUG) << "@@Flux = " << flx_driver->PdgCode() << " , E =" << flx_driver->Momentum().E() << " GeV";
    LOG("ComputeAttenuation", pDEBUG) << "Postion   = [ " << flx_driver->Position().X() << " m, " << flx_driver->Position().Y() << " m, " << flx_driver->Position().Z() << " m, " << flx_driver->Position().T() << " s ] ";
    LOG("ComputeAttenuation", pDEBUG) << "Direction = [ " << fdir.X() << " , " << fdir.Y() << " , " << fdir.Z() << " ]";

    int eventid = trj_driver->GenerateEvent();

    if ( eventid==1 ) {

      EventRecord * event = trj_driver->GetEvent();

      const TLorentzVector & vtx  = *event->Vertex();
      TVector3 nudir = event->Probe()->P4()->Vect().Unit();
 
      // print-out
      LOG("ComputeAttenuation", pDEBUG) << "Neutrino interaction!!!";
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << vtx.X() << " m, " << vtx.Y() << " m, " << vtx.Z() << " m, " << vtx.T() << " s ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << nudir.X() << " , " << nudir.Y() << " , " << nudir.Z() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << *event;
      assert(sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y()+vtx.Z()*vtx.Z())<gOptGeoLimit);

      NInt++;

      GHepParticle * p = 0;
      TObjArrayIter piter(event);

      while ( (p=(GHepParticle *)piter.Next()) ) {
        
        //quit when it reachs the hadronic shower to avoid counting neutrinos from top
        //if ( abs(p->Pdg())==2000000001 ) break; 

        if ( p->E()<=spline_Erange.min || p->Status()!=kIStStableFinalState ) continue;

        p->SetPosition( vtx.X()+p->Vx()*1e-15, vtx.Y()+p->Vy()*1e-15, vtx.Z()+p->Vz()*1e-15, vtx.T()+p->Vt() ); //position -> passing from fm(genie) to m(geom)

        if ( pdg::IsNeutrino(TMath::Abs(p->Pdg())) ) {
          p->SetEnergy( TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()) ); //not using E to avoid problems with neutrinos from pythia not on shell
          SecNu.push_back(*p); // X -> nu
        }

        else if ( pdg::IsTau(TMath::Abs(p->Pdg())) ) {
          std::vector<GHepParticle> TauProd = tauprop->Propagate(p,spline_Erange.min);
          for ( auto & tpr : TauProd ) {
            if ( tpr.E()>spline_Erange.min ) {                
              if      ( pdg::IsTau     (TMath::Abs(tpr.Pdg())) ) OutPart.push_back(tpr); // X -> tau
              else if ( pdg::IsNeutrino(TMath::Abs(tpr.Pdg())) ) SecNu.push_back(tpr);   // X -> tau -> nu
            }
          }
        }

        else if (gOptHadronProp) {

          TParticlePDG * partinfo = PdgDB->GetParticle(p->Pdg());
          string pclass = partinfo->ParticleClass();

          if( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) { 
            std::vector<GHepParticle> CharmProd = hadronprop->Propagate(p,pclass);
            for ( auto & cpr : CharmProd ) {
              if ( cpr.E()>spline_Erange.min ) {  
                if ( pdg::IsTau(TMath::Abs(cpr.Pdg())) ) {
                  std::vector<GHepParticle> TauProd = tauprop->Propagate(&cpr,spline_Erange.min);
                  for ( auto & tpr : TauProd ) {
                    if ( tpr.E()>spline_Erange.min ) {  
                      if      ( pdg::IsTau     (TMath::Abs(tpr.Pdg())) ) OutPart.push_back(tpr); // X -> charmed -> tau
                      else if ( pdg::IsNeutrino(TMath::Abs(tpr.Pdg())) ) SecNu.push_back(tpr);   // X -> charmed -> tau -> nu
                    }
                  }
                }
                else if ( pdg::IsNeutrino(TMath::Abs(cpr.Pdg())) ) SecNu.push_back(cpr); // X -> charmed -> nu
                else OutPart.push_back(cpr); // X -> charmed
              }
            }
          }
          
          else if( pclass=="B-Meson" || pclass=="B-Baryon" ) { 
            
            std::vector<GHepParticle> BottomProd = hadronprop->Propagate(p,pclass);
            for ( auto & bpr : BottomProd ) {
              if ( bpr.E()>spline_Erange.min ) {  
                if ( pdg::IsTau(TMath::Abs(bpr.Pdg())) ) {
                  std::vector<GHepParticle> TauProd = tauprop->Propagate(&bpr,spline_Erange.min);
                  for ( auto & tpr : TauProd ) {
                    if ( tpr.E()>spline_Erange.min ) {  
                      if      ( pdg::IsTau     (TMath::Abs(tpr.Pdg())) ) OutPart.push_back(tpr); // X -> bottom -> tau
                      else if ( pdg::IsNeutrino(TMath::Abs(tpr.Pdg())) ) SecNu.push_back(tpr);   // X -> bottom -> tau -> nu
                    }
                  }
                }
                else if ( pdg::IsNeutrino(TMath::Abs(bpr.Pdg())) ) SecNu.push_back(bpr); // X -> bottom -> nu
                else {
                  partinfo = PdgDB->GetParticle(bpr.Pdg());
                  pclass = partinfo->ParticleClass();
                  if( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) { 
                    std::vector<GHepParticle> CharmProd = hadronprop->Propagate(&bpr,pclass);
                    for ( auto & cpr : CharmProd ) {
                      if ( cpr.E()>spline_Erange.min ) {  
                        if ( pdg::IsTau(TMath::Abs(cpr.Pdg())) ) {
                          std::vector<GHepParticle> TauProd = tauprop->Propagate(&cpr,spline_Erange.min);
                          for ( auto & tpr : TauProd ) {
                            if ( tpr.E()>spline_Erange.min ) {  
                              if      ( pdg::IsTau     (TMath::Abs(tpr.Pdg())) ) OutPart.push_back(tpr); // X -> bottom -> charm -> tau
                              else if ( pdg::IsNeutrino(TMath::Abs(tpr.Pdg())) ) SecNu.push_back(tpr);   // X -> bottom -> charm -> tau -> nu
                            }
                          }
                        }
                        else if ( pdg::IsNeutrino(TMath::Abs(cpr.Pdg())) ) SecNu.push_back(cpr); // X -> bottom -> charm  -> nu
                        else OutPart.push_back(cpr); // X -> bottom -> charm 
                      }
                    }
                  }
                  else OutPart.push_back(bpr); // X -> bottom
                }
              }
            }
          } //bottom prop
        
        } //hadron prop


      }

      delete event;

    }
    else if ( eventid==2 ) {
      LOG("ComputeAttenuation", pDEBUG) << "Neutrino exiting!!!";
      OutPart.push_back(*flx_driver->GetNeutrino());      
    }
    else {
      LOG("ComputeAttenuation", pWARN) << "Failed trial: " << trial;
      if ( trial==NMAXTRIALS-1 ){
        LOG("ComputeAttenuation", pFATAL) << "Neutrinos can not be be generated" ;
        LOG("ComputeAttenuation", pFATAL) << "Trials: " << trial ;
        exit(1);
      }      
      SecNu.clear();
      OutPart.clear();
      NInt = 0;
      trial++;
      continue;
    }

    if ( SecNu.size()==0 ) {
      LOG("ComputeAttenuation", pDEBUG) << "No more secondary neutrinos!!!";
      influx->Fill();
      LOG("ComputeAttenuation", pINFO) << "Summary: ";
      LOG("ComputeAttenuation", pINFO) << In_Pdg << " ----> ((( )))";
      for (unsigned int i=0; i<OutPart.size(); i++) {
        FillParticle(&OutPart[i],Out_Pdg,Out_PX,Out_PY,Out_PZ,Out_E,Out_X,Out_Y,Out_Z);
        outflux->Fill();
        LOG("ComputeAttenuation", pINFO) << "         ((( ))) ----> " << Out_Pdg;
      }
      LOG("ComputeAttenuation", pINFO) << "Total In: " << influx->GetEntries();
      LOG("ComputeAttenuation", pINFO) << "Total Out: " << outflux->GetEntries();
      OutPart.clear();
      NInt = 0;
      trial = 0;
      NNu++;
    }

  } 

  influx->Write();
  outflux->Write();
  outfile->Close();

  // clean-up
  delete flx_driver;
  delete trj_driver;

  return 0;
}

//**************************************************************************
//**************************************************************************
//THIS FUNCTION READS THE INPUT OPTIONS
//**************************************************************************
//**************************************************************************
void GetCommandLineArgs(int argc, char ** argv)
{

  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  size_t apos;
  string opt;

  for (int i = 1; i < argc; i++) {
    string arg(argv[i]);
    if (arg.compare(0,1,"-")==0){
      apos=arg.find(" ");
      opt=arg.substr(0,apos);

      if(opt.compare("--number-of-events")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading number of events to generate";
        gOptNev = atof(argv[i]);
      }          
      if(opt.compare("--geometry")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading geometry";
        gOptGeometry = argv[i];
      }          
      if(opt.compare("--geometry-limit")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading geometry limit";
        gOptGeoLimit = atof(argv[i]);
      }          
      if(opt.compare("--seed")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading seed";
        gOptRanSeed = atoi(argv[i]);
      }          
      if(opt.compare("--output")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading output name";
        gOptOutName = argv[i];
      }          
      if(opt.compare("--probe")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino flavor";
        gOptPdg = atoi(argv[i]);
      }          
      if(opt.compare("--alpha")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino spectrum alpha";
        gOptAlpha = atof(argv[i]);
      }          
      if(opt.compare("--costheta")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino angle";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if      (saux.size()==1) {
          gOptCthmin = gOptCthmax = atof(saux[0].c_str());
        }
        else if (saux.size()==2) {
          gOptCthmin = atof(saux[0].c_str());
          gOptCthmax = atof(saux[1].c_str());
        }
        else {
          LOG("ComputeAttenuation", pFATAL) << "Wrong costheta argument " << argv[i] <<"; exit";
          exit(1);
        }
      }          
      if(opt.compare("--energy")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino min energy";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if      (saux.size()==1) {
          gOptEmin = gOptEmax = atof(saux[0].c_str());
        }
        else if (saux.size()==2) {
          gOptEmin = atof(saux[0].c_str());
          gOptEmax = atof(saux[1].c_str());
        }
        else {
          LOG("ComputeAttenuation", pFATAL) << "Wrong energy argument " << argv[i] << "; exit";
          exit(1);
        }
      }          
      if(opt.compare("--cross-sections")==0){ 
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading cross-section file";
        gOptInpXSecFile = argv[i];  
      }
      if(opt.compare("--detector-position")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector pos";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if (saux.size()==3) {
          gOptDetPos[0] = atof(saux[0].c_str());
          gOptDetPos[1] = atof(saux[1].c_str());
          gOptDetPos[2] = atof(saux[2].c_str());
        }
        else {
          LOG("ComputeAttenuation", pFATAL) << "Wrong detector argument " << argv[i] << "; exit";
          exit(1);
        }
      }       
      if(opt.compare("--detector-radius")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector radius";
        gOptDetRadius = atof(argv[i]);  
      }
      if(opt.compare("--detector-height")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector height";
        gOptDetHeight = atof(argv[i]);  
      }       
      if(opt.compare("--offset")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading offset";
        gOptOffset = atof(argv[i]);  
      }       
      if(opt.compare("--tau-propagation")==0){ 
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading tau propagation tables";
        gOptTauProp = argv[i];  
      }      
      if(opt.compare("--no-hadron-propagation")==0){ 
        LOG("ComputeAttenuation", pINFO) << "Disable hadron propagation";
        gOptHadronProp = false;  
      }      

    }
  }



  if ( gSystem->AccessPathName(gOptGeometry.c_str()) ) {
    LOG("ComputeAttenuation", pFATAL) << "Geometry file not found: " << gOptGeometry << "; exit";
    exit(1);
  }

  if ( gSystem->AccessPathName(gOptInpXSecFile.c_str()) ) {
    LOG("ComputeAttenuation", pFATAL) << "Cross section file not found: " << gOptInpXSecFile << "; exit";
    exit(1);
  }

  if ( !pdg::IsNeutrino(abs(gOptPdg)) && gOptPdg!=-1 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong neutrino flavor: " << gOptPdg << "; exit";
    exit(1);
  }

  if ( gOptNev<=0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong number of events: " << gOptNev << "; exit";
    exit(1);    
  }

  if ( gOptGeoLimit<=0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong geomtry limit: " << gOptGeoLimit << "; exit";
    exit(1);    
  }

  if ( gOptDetHeight<0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong height: " << gOptDetHeight << "; exit";
    exit(1);    
  }

  if ( gOptDetRadius<0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong radius: " << gOptDetRadius << "; exit";
    exit(1);    
  }

  if ( gOptOffset<0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong offset: " << gOptOffset << "; exit";
    exit(1);    
  }

  if (gOptEmin>gOptEmax || gOptEmin<spline_Erange.min || gOptEmax>spline_Erange.max) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong energy range: " << gOptEmin << "," << gOptEmax << "; exit";
    exit(1);    
  }
  if (gOptEmin==gOptEmax) LOG("ComputeAttenuation", pNOTICE) << "Running in monoenergetic mode";

  if (gOptCthmin>gOptCthmax || gOptCthmin<-1. || gOptCthmax<-1. || gOptCthmin>1. || gOptCthmax>1. ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong costheta range: " << gOptCthmin << "," << gOptCthmax << "; exit";
    exit(1);    
  }
  if (gOptCthmin==gOptCthmax) LOG("ComputeAttenuation", pNOTICE) << "Running in monoangular mode";

  if ( gSystem->AccessPathName(gOptTauProp.c_str()) ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong tau propagation: " << gOptTauProp << "; exit";
    exit(1);          
  }


  LOG("ComputeAttenuation", pNOTICE) << "\n\n" << utils::print::PrintFramedMesg("NuPropEarth job configuration");
  LOG("ComputeAttenuation", pNOTICE) 
     << "\n"
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using cross-section file: " << gOptInpXSecFile
     << "\n @@ Using output file name: " << gOptOutName
     << "\n @@ Geometry"
     << "\n\t File          = " << gOptGeometry
     << "\n\t DetPos        = " << gOptDetPos[0] << "," << gOptDetPos[1] << "," << gOptDetPos[2] << " m"
     << "\n\t DetRadius     = " << gOptDetRadius << " m"
     << "\n\t DetHeight     = " << gOptDetHeight << " m"
     << "\n\t Offset        = " << gOptOffset << " m"
     << "\n @@ Exposure" 
     << "\n\t" << gOptNev
     << "\n @@ Kinematics" 
     << "\n\t Pdg          = " << gOptPdg
     << "\n\t CosTheta_min = " << gOptCthmin
     << "\n\t CosTheta_max = " << gOptCthmax
     << "\n\t Emin         = " << gOptEmin << " GeV"
     << "\n\t Emax         = " << gOptEmax << " GeV"
     << "\n\t Alpha        = " << gOptAlpha
     << "\n @@ Propagation" 
     << "\n\t TauProp      = " << gOptTauProp
     << "\n\t HadronProp   = " << gOptHadronProp
     << "\n\n";


}



