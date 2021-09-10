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
void FillParticle       (GHepParticle * part, int &pdg, double &px, double &py,double &pz, double &e, double &vx, double &vy, double &vz, double &t) {
  pdg = part->Pdg();
  px  = part->Px(); py = part->Py(); pz = part->Pz(); e = part->E();
  vx  = part->Vx(); vy = part->Vy(); vz = part->Vz(); t = part->Vt();
}


// User-specified options:
string          gOptOutName = "./test.root";               
string          gOptGeometry = "";               
int             gOptRanSeed = 0;               
string          gOptInpXSecFile = "";           
int             gOptNev = 0;              
int             gOptPdg = -1; //all flavors                   
double          gOptEmin = 1e2;
double          gOptEmax = 1e10;
double          gOptCthmin = -1.;
double          gOptCthmax =  1.;
double          gOptAlpha = 1; //flat spectrum in log10e
double          gOptDetPos[3] = { 0., 0., 0. }; //Detector position at center of the volume
double          gOptRadius = 0.; //Detector Radius (radius of a cylindrical shape) originally a point
double          gOptHeight = 0.; //Detector height (height of the cylinder centerd at gOptDetPos[2]) originally a point
string          gOptTauProp = ""; //"","NOELOSS","TAUSIC-ALLM","TAUSIC-BS","PROPOSAL"


const int NMAXINT = 100;
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

  GeomAnalyzerI * geom_driver = dynamic_cast<GeomAnalyzerI *> (rgeom);

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Tau Propagation...";
  TauPropagation * tauprop = new TauPropagation(gOptTauProp,gOptRanSeed,geom_driver);

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Hadron Propagation...";
  HadronPropagation * hadronprop = new HadronPropagation(geom_driver);
  TDatabasePDG * PdgDB  = TDatabasePDG::Instance();



  LOG("ComputeAttenuation", pDEBUG) << "Creating GFluxI...";
  IncomingFlux * flx_driver = new IncomingFlux(gOptPdg, gOptAlpha, gOptCthmin, gOptCthmax, gOptEmin, gOptEmax, gOptDetPos, gOptRadius, gOptHeight);

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
  trj_driver->UseGeomAnalyzer(geom_driver);

  LOG("ComputeAttenuation", pDEBUG) << "Configuring GTRJDriver...";
  trj_driver->Configure(spline_Erange.min,spline_Erange.max);

  // create file so tree is saved inside
  TFile * outfile = new TFile(gOptOutName.c_str(),"RECREATE");
  
  int In_Pdg,Out_Pdg;
  double In_X4[4],Out_X4[4];
  double In_P4[4],Out_P4[4];
  int NInt,NIntCC,NIntNC;
  int SctID[NMAXINT];
  int IntID[NMAXINT];
  int Tgt[NMAXINT];
  int Pdg[NMAXINT];

  TTree *influx = new TTree("InFlux","InFlux");
  influx->Branch("In_Pdg", &In_Pdg, "In_Pdg/I"   );       
  influx->Branch("In_X4",  In_X4,   "In_X4[4]/D" );       
  influx->Branch("In_P4",  In_P4,   "In_P4[4]/D" );       

  TTree *outflux = new TTree("OutFlux","OutFlux");
  outflux->Branch("In_Pdg",  &In_Pdg,  "In_Pdg/I"      );       
  outflux->Branch("In_X4",   In_X4,    "In_X4[4]/D"    );       
  outflux->Branch("In_P4",   In_P4,    "In_P4[4]/D"    );       
  outflux->Branch("Out_Pdg", &Out_Pdg, "Out_Pdg/I"     );       
  outflux->Branch("Out_X4",  Out_X4,   "Out_X4[4]/D"   );       
  outflux->Branch("Out_P4",  Out_P4,   "Out_P4[4]/D"   );       
  outflux->Branch("NInt",    &NInt,    "NInt/I"        );       
  outflux->Branch("NIntCC",  &NIntCC,  "NIntCC/I"      );       
  outflux->Branch("NIntNC",  &NIntNC,  "NIntNC/I"      );       
  outflux->Branch("SctID",   SctID,    "SctID[NInt]/I" );       
  outflux->Branch("IntID",   IntID,    "IntID[NInt]/I" );       
  outflux->Branch("Tgt",     Tgt,      "Tgt[NInt]/I"   );        
  outflux->Branch("Pdg",     Pdg,      "Pdg[NInt]/I"   );       

  LOG("ComputeAttenuation", pDEBUG) << "Event loop...";

  GHepParticle * Nu_In = new GHepParticle();
  GHepParticle * Nu_Out = new GHepParticle();
  std::vector<GHepParticle> SecNu;
  std::vector<GHepParticle> OutTau;

  int NNu = 0;
  NIntCC = 0;
  NIntNC = 0;
  NInt = 0;
  while ( NNu<gOptNev ) {

    if ( SecNu.size()>0 ) {
      LOG("ComputeAttenuation", pDEBUG) << "@@SecNu      = "   << SecNu[0].Pdg() << " , E = " << SecNu[0].E() << " GeV";
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << SecNu[0].Vx() << " m, " << SecNu[0].Vy() << " m, " << SecNu[0].Vz() << " m, " << SecNu[0].Vt() << " s ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << SecNu[0].Px()/SecNu[0].E() << " , " << SecNu[0].Py()/SecNu[0].E() << " , " << SecNu[0].Pz()/SecNu[0].E() << " ]";
      flx_driver->InitNeutrino(SecNu[0]);
      SecNu.erase(SecNu.begin()); //remove first entry because it is just process
    }

    EventRecord * event = trj_driver->GenerateEvent();

    if (NInt==0) {
      if ( NNu%100==0 ) LOG("ComputeAttenuation", pINFO) << "Event " << NNu << " out of " << gOptNev;
      Nu_In = flx_driver->GetNeutrino();
      FillParticle(Nu_In,In_Pdg,In_P4[0],In_P4[1],In_P4[2],In_P4[3],In_X4[0],In_X4[1],In_X4[2],In_X4[3]);
      influx->Fill();
    }

    LOG("ComputeAttenuation", pDEBUG) << "Shooting Neutrino: " << NNu << "( " << NInt << " ) ----->" ;
    LOG("ComputeAttenuation", pDEBUG) << "@@Flux = " << flx_driver->GetNeutrino()->Pdg() << " , E =" << flx_driver->GetNeutrino()->E() << " GeV";
    LOG("ComputeAttenuation", pDEBUG) << "Postion   = [ " << flx_driver->GetNeutrino()->Vx() << " m, " << flx_driver->GetNeutrino()->Vy() << " m, " << flx_driver->GetNeutrino()->Vz() << " m, " << flx_driver->GetNeutrino()->Vt() << " s ] ";
    LOG("ComputeAttenuation", pDEBUG) << "Direction = [ " << flx_driver->GetNeutrino()->Px()/flx_driver->GetNeutrino()->E() << " , " << flx_driver->GetNeutrino()->Py()/flx_driver->GetNeutrino()->E() << " , " << flx_driver->GetNeutrino()->Pz()/flx_driver->GetNeutrino()->E() << " ]";

    if (event) {

      SctID[NInt] = event->Summary()->ProcInfoPtr()->ScatteringTypeId();
      IntID[NInt] = event->Summary()->ProcInfoPtr()->InteractionTypeId();
      Tgt[NInt]   = ( event->TargetNucleus() ) ? event->TargetNucleus()->Pdg() : event->HitNucleon()->Pdg();
      Pdg[NInt]   = event->Probe()->Pdg();

      const TLorentzVector & X4  = *event->Vertex();
      const TLorentzVector & P4  = *event->Probe()->P4();

      // print-out
      LOG("ComputeAttenuation", pDEBUG) << "Neutrino interaction!!!";
      LOG("ComputeAttenuation", pDEBUG) << "@@Interact   = "   << SctID[NInt] << "   " << IntID[NInt];
      LOG("ComputeAttenuation", pDEBUG) << "@@Target     = "   << Tgt[NInt];
      LOG("ComputeAttenuation", pDEBUG) << "@@Probe      = "   << Pdg[NInt] << " , E =" << P4.E() << " GeV";
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << X4.X() << " m, " << X4.Y() << " m, " << X4.Z() << " m, " << X4.T() << " s ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << P4.Px()/P4.E() << " , " << P4.Py()/P4.E() << " , " << P4.Pz()/P4.E() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << *event;
      
      if (IntID[NInt]==2) NIntCC++;
      else                NIntNC++;
      NInt++;

      GHepParticle * p = 0;
      TObjArrayIter piter(event);

      while ( (p=(GHepParticle *)piter.Next()) ) {
        
        //quit when it reachs the hadronic shower to avoid counting neutrinos from top
        //if ( abs(p->Pdg())==2000000001 ) break; 

        if ( p->E()<=spline_Erange.min || p->Status()!=kIStStableFinalState ) continue;

        if ( pdg::IsNeutrino(TMath::Abs(p->Pdg())) ) {
          LOG("ComputeAttenuation", pDEBUG) << p->Pdg() << " " << p->E();
          p->SetEnergy( TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()) ); //not using E to avoid problems with neutrinos from pythia not on shell
          p->SetPosition( X4.X()+p->Vx()*1e-15, X4.Y()+p->Vy()*1e-15, X4.Z()+p->Vz()*1e-15, X4.T()+p->Vt() ); //position -> passing from fm(genie) to m(geom)
          SecNu.push_back(*p);
        }
        else if ( pdg::IsTau(TMath::Abs(p->Pdg())) ) {
          p->SetPosition( X4.X()+p->Vx()*1e-15, X4.Y()+p->Vy()*1e-15, X4.Z()+p->Vz()*1e-15, X4.T()+p->Vt() ); //position -> passing from fm(genie) to m(geom)
          std::vector<GHepParticle> TauProd = tauprop->Propagate(p);
          for ( auto & tpr : TauProd ) {
            if ( tpr.E()>spline_Erange.min ) {  
              if ( pdg::IsTau(TMath::Abs(tpr.Pdg())) ) OutTau.push_back(tpr);
              else SecNu.push_back(tpr);
            }
          }
        }
        else {
          TParticlePDG * partinfo = PdgDB->GetParticle(p->Pdg());
          string pclass = partinfo->ParticleClass();
          if( pclass=="B-Meson" || pclass=="B-Baryon" || pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) { 
            p->SetPosition( X4.X()+p->Vx()*1e-15, X4.Y()+p->Vy()*1e-15, X4.Z()+p->Vz()*1e-15, X4.T()+p->Vt() ); //position -> passing from fm(genie) to m(geom)
            std::vector<GHepParticle> HadronProd = hadronprop->Propagate(p,pclass);
            for ( auto & hpr : HadronProd ) {
              if ( hpr.E()>spline_Erange.min ) {  
                if ( pdg::IsTau(TMath::Abs(hpr.Pdg())) ) {
                  cout << "One tau from hadron: " << p->Pdg() << endl;
                  std::vector<GHepParticle> TauProd = tauprop->Propagate(&hpr);
                  for ( auto & tpr : TauProd ) {
                    if ( tpr.E()>spline_Erange.min ) {  
                      if ( pdg::IsTau(TMath::Abs(tpr.Pdg())) ) OutTau.push_back(tpr);
                      else SecNu.push_back(tpr);
                    }
                  }
                }
                else SecNu.push_back(hpr);
              }
            }
          }
        }

      }

      delete event;

    }
    else {
      LOG("ComputeAttenuation", pDEBUG) << " ----> Neutrino: Goodbye Earth!!!";
      Nu_Out = flx_driver->GetNeutrino();
      FillParticle(Nu_Out,Out_Pdg,Out_P4[0],Out_P4[1],Out_P4[2],Out_P4[3],Out_X4[0],Out_X4[1],Out_X4[2],Out_X4[3]);
      outflux->Fill();
    }

    if ( SecNu.size()==0 ) {
      LOG("ComputeAttenuation", pDEBUG) << " ----> No more secondary neutrinos";
      for (unsigned int i=0; i<OutTau.size(); i++) {
        LOG("ComputeAttenuation", pDEBUG) << " ----> Tau: Goodbye Earth!!!";
        FillParticle(&OutTau[i],Out_Pdg,Out_P4[0],Out_P4[1],Out_P4[2],Out_P4[3],Out_X4[0],Out_X4[1],Out_X4[2],Out_X4[3]);
        outflux->Fill();
      }
      OutTau.clear();
      NInt = 0;
      NIntCC = 0;
      NIntNC = 0;
      NNu++;
    }

  } 

  influx->Write();
  outflux->Write();
  outfile->Close();

  // clean-up
  delete geom_driver;
  delete flx_driver;
  delete trj_driver;

  return 0;
}

//**************************************************************************
//**************************************************************************
//THIS FUNCTION DECAY TAUS
//**************************************************************************
//**************************************************************************


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
      if(opt.compare("--radius")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector radius";
        gOptRadius = atof(argv[i]);  
      }
      if(opt.compare("--height")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector height";
        gOptHeight = atof(argv[i]);  
      }       
      if(opt.compare("--tau-propagation")==0){ 
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Enable tau propagation";
        gOptTauProp = argv[i];  
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

  if ( gOptHeight<0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong height: " << gOptHeight << "; exit";
    exit(1);    
  }

  if ( gOptRadius<0 ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong radius: " << gOptRadius << "; exit";
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

  if (gOptTauProp!="" && gOptTauProp!="NOELOSS" && gOptTauProp!="TAUSIC-ALLM" && gOptTauProp!="TAUSIC-BS" && gOptTauProp!="PROPOSAL" ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong tau propagation: " << gOptTauProp << "; exit";
    exit(1);    
  }


  LOG("ComputeAttenuation", pNOTICE) << "\n\n" << utils::print::PrintFramedMesg("NuPropEarth job configuration");
  LOG("ComputeAttenuation", pNOTICE) 
     << "\n"
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using cross-section file: " << gOptInpXSecFile
     << "\n @@ Using output file name: " << gOptOutName
     << "\n @@ Exposure" 
     << "\n\t" << gOptNev
     << "\n @@ Kinematics" 
     << "\n\t Pdg          = " << gOptPdg
     << "\n\t CosTheta_min = " << gOptCthmin
     << "\n\t CosTheta_max = " << gOptCthmax
     << "\n\t Emin         = " << gOptEmin << " GeV"
     << "\n\t Emax         = " << gOptEmax << " GeV"
     << "\n\t Alpha        = " << gOptAlpha
     << "\n\t DetPos       = " << gOptDetPos[0] << "," << gOptDetPos[1] << "," << gOptDetPos[2] << " m"
     << "\n\t Radius       = " << gOptRadius << " m"
     << "\n\t Height       = " << gOptHeight << " m"
     << "\n\t TauProp      = " << gOptTauProp
     << "\n\n";


}



