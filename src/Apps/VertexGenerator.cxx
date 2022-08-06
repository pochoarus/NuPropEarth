#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TNamed.h>

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


using std::string;
using std::pair;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

//functions
void   GetCommandLineArgs (int argc, char ** argv);
double DistToBorder       (double vx, double vy, double vz, double dx, double dy, double dz, double radius, double lw, double up);

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

const Range1D_t spline_Erange = { 1e2, 1e12 }; //energy range limit based on HEDIS splines 

//**************************************************************************
//**************************************************************************
//THIS IS THE MAIN FUNCTION
//**************************************************************************
//**************************************************************************
int main(int argc, char** argv)
{

  // Parse command line arguments
  LOG("VertexGenerator", pDEBUG) << "Reading options...";
  GetCommandLineArgs(argc,argv);

  geometry::ROOTGeomAnalyzer * rgeom = new geometry::ROOTGeomAnalyzer(gOptGeometry);
  rgeom -> SetLengthUnits  ( genie::utils::units::UnitFromString("m")     );
  rgeom -> SetDensityUnits ( genie::utils::units::UnitFromString("g_cm3") );
  rgeom -> SetTopVolName   ("");

  TGeoVolume * topvol = rgeom->GetGeometry()->GetTopVolume();
   if(!topvol) {
    LOG("VertexGenerator", pFATAL) << " ** Null top ROOT geometry volume!";
    exit(1);
  }
  LOG("VertexGenerator", pNOTICE) << topvol->GetName();

  TDatabasePDG * PdgDB  = TDatabasePDG::Instance();

  LOG("VertexGenerator", pDEBUG) << "Creating GFluxI...";
  IncomingFlux * flx_driver = new IncomingFlux(gOptPdg, gOptAlpha, gOptCthmin, gOptCthmax, gOptEmin, gOptEmax, gOptDetPos, gOptDetRadius, gOptDetHeight, gOptOffset, true);

  RunOpt::Instance()->BuildTune();

  // Iinitialization of random number generators, cross-section table, messenger, cache etc...
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  LOG("VertexGenerator", pDEBUG) << "Creating GTRJDriver...";
  GTRJDriver * trj_driver = new GTRJDriver;
  trj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  trj_driver->UseFluxDriver(dynamic_cast<GFluxI *>(flx_driver));
  trj_driver->UseGeomAnalyzer(rgeom);

  LOG("VertexGenerator", pDEBUG) << "Configuring GTRJDriver...";
  trj_driver->Configure(spline_Erange.min,spline_Erange.max);

  // create file so tree is saved inside
  TFile * outfile = new TFile(gOptOutName.c_str(),"RECREATE");
  
  double Gun_X,Gun_Y,Gun_Z;
  double Wgt;
  vector<int>    SctID;
  vector<int>    IntID;
  vector<int>    Tgt;
  vector<double> X,Y,Z,T;
  vector<int>    In_Pdg;
  vector<double> In_PX,In_PY,In_PZ,In_E;
  vector<vector<double>> Out_Pdg;
  vector<vector<double>> Out_PX,Out_PY,Out_PZ,Out_E;

  TTree *outflux = new TTree("OutFlux","OutFlux");
  outflux->Branch("Gun_X",    &Gun_X,  "Gun_X/D" );       
  outflux->Branch("Gun_Y",    &Gun_Y,  "Gun_Y/D" );       
  outflux->Branch("Gun_Z",    &Gun_Z,  "Gun_Z/D" );       
  outflux->Branch("Wgt",      &Wgt,   "Wgt/D"    );       
  outflux->Branch("SctID",    &SctID   );
  outflux->Branch("IntID",    &IntID   );
  outflux->Branch("Tgt",      &Tgt     );
  outflux->Branch("X",        &X       );
  outflux->Branch("Y",        &Y       );
  outflux->Branch("Z",        &Z       );
  outflux->Branch("T",        &T       );
  outflux->Branch("In_Pdg",   &In_Pdg  );
  outflux->Branch("In_PX",    &In_PX   );
  outflux->Branch("In_PY",    &In_PY   );
  outflux->Branch("In_PZ",    &In_PZ   );
  outflux->Branch("In_E",     &In_E    );
  outflux->Branch("Out_Pdg",  &Out_Pdg );
  outflux->Branch("Out_PX",   &Out_PX  );
  outflux->Branch("Out_PY",   &Out_PY  );
  outflux->Branch("Out_PZ",   &Out_PZ  );
  outflux->Branch("Out_E",    &Out_E   );

  LOG("VertexGenerator", pDEBUG) << "Event loop...";

  std::vector<GHepParticle> SecNu;

  int NNu  = 0;
  int NInt = 0;
  while ( NNu<gOptNev ) {

    bool force = false;
    if ( SecNu.size()==0 ) {
      if ( NNu%100==0 ) LOG("VertexGenerator", pNOTICE) << "Event " << NNu << " out of " << gOptNev;
      flx_driver->InitNeutrino();
      Gun_X = flx_driver->Position().X();
      Gun_Y = flx_driver->Position().Y();
      Gun_Z = flx_driver->Position().Z();
      Wgt   = flx_driver->GetWeight(flx_driver->Momentum().E());
      force = true;
    }
    else {
      LOG("VertexGenerator", pDEBUG) << "@@SecNu = "   << SecNu[0].Pdg() << " , E = " << SecNu[0].E() << " GeV";
      flx_driver->InitNeutrino(SecNu[0]);
      SecNu.erase(SecNu.begin()); //remove first entry because it is just process
    }

    TVector3 fdir = flx_driver->Momentum().Vect().Unit();
    double detx = flx_driver->Position().X()-gOptDetPos[0];
    double dety = flx_driver->Position().Y()-gOptDetPos[1];
    double detz = flx_driver->Position().Z()-gOptDetPos[2];
    double dist = DistToBorder(detx,dety,detz,fdir.X(),fdir.Y(),fdir.Z(),gOptDetRadius,-gOptDetHeight/2.,gOptDetHeight/2.);
    assert(dist>=0.);

    LOG("VertexGenerator", pDEBUG) << "Shooting Neutrino: " << NNu << "( " << NInt << " ) ----->" ;
    LOG("VertexGenerator", pDEBUG) << "@@Flux = " << flx_driver->PdgCode() << " , E =" << flx_driver->Momentum().E() << " GeV";
    LOG("VertexGenerator", pDEBUG) << "Postion   = [ " << flx_driver->Position().X()-gOptDetPos[0] << " m, " << flx_driver->Position().Y()-gOptDetPos[1] << " m, " << flx_driver->Position().Z()-gOptDetPos[2] << " m, " << flx_driver->Position().T() << " s ] ";
    LOG("VertexGenerator", pDEBUG) << "Direction = [ " << fdir.X() << " , " << fdir.Y() << " , " << fdir.Z() << " ]";
    LOG("VertexGenerator", pDEBUG) << "dist      = " << dist;

    int eventid = trj_driver->GenerateEvent(force,dist);

    if ( eventid==1 ) {

      EventRecord * event = trj_driver->GetEvent();
      Wgt *= trj_driver->GetWeight();

      const TLorentzVector & vtx  = *event->Vertex();
      TVector3 nudir = event->Probe()->P4()->Vect().Unit();

      SctID .push_back(event->Summary()->ProcInfoPtr()->ScatteringTypeId());
      IntID .push_back(event->Summary()->ProcInfoPtr()->InteractionTypeId());
      Tgt   .push_back(( event->TargetNucleus() ) ? event->TargetNucleus()->Pdg() : event->HitNucleon()->Pdg());      
      In_Pdg.push_back(event->Probe()->Pdg());
      In_PX .push_back(event->Probe()->Px());
      In_PY .push_back(event->Probe()->Py());
      In_PZ .push_back(event->Probe()->Pz());
      In_E  .push_back(event->Probe()->E());
      X     .push_back(vtx.X());
      Y     .push_back(vtx.Y());
      Z     .push_back(vtx.Z());
      T     .push_back(vtx.T());

      // print-out
      LOG("VertexGenerator", pDEBUG) << "Neutrino interaction!!!";
      LOG("VertexGenerator", pDEBUG) << SctID[NInt] << " " << IntID[NInt];
      LOG("VertexGenerator", pDEBUG) << "  Position   = [ " << vtx.X()-gOptDetPos[0] << " m, " << vtx.Y()-gOptDetPos[1] << " m, " << vtx.Z()-gOptDetPos[2] << " m, " << vtx.T() << " s ]";
      LOG("VertexGenerator", pDEBUG) << "  Direction  = [ " << nudir.X() << " , " << nudir.Y() << " , " << nudir.Z() << " ]";
      LOG("VertexGenerator", pDEBUG) << *event;
      assert(sqrt(vtx.X()*vtx.X()+vtx.Y()*vtx.Y()+vtx.Z()*vtx.Z())<gOptGeoLimit);

      vector<double> Aux_Pdg;
      vector<double> Aux_PX,Aux_PY,Aux_PZ,Aux_E;

      GHepParticle * p = 0;
      TObjArrayIter piter(event);
      while ( (p=(GHepParticle *)piter.Next()) ) {
        
        if ( p->Status()!=kIStStableFinalState ) continue;

        p->SetPosition( vtx.X()+p->Vx()*1e-15, vtx.Y()+p->Vy()*1e-15, vtx.Z()+p->Vz()*1e-15, vtx.T()+p->Vt() ); //position -> passing from fm(genie) to m(geom)

        bool outfill = false;

        if ( pdg::IsNeutrino(TMath::Abs(p->Pdg())) ) {
          outfill = true;
          if ( p->E()>spline_Erange.min ) {
            p->SetEnergy( TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()) ); //not using E to avoid problems with neutrinos from pythia not on shell
            SecNu.push_back(*p); // X -> nu
          }
        }
        else if ( pdg::IsChargedLepton(TMath::Abs(p->Pdg())) ) outfill = true;
        // else {
        //   TParticlePDG * partdb = PdgDB->GetParticle(p->Pdg());
        //   string pclass = partdb->ParticleClass();
        //   if      ( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) outfill = true;
        //   else if ( pclass==     "B-Meson" || pclass==     "B-Baryon" ) outfill = true;
        // }

        if (outfill) {
          Aux_Pdg.push_back(p->Pdg());
          Aux_PX.push_back(p->Px());
          Aux_PY.push_back(p->Py());
          Aux_PZ.push_back(p->Pz());
          Aux_E.push_back(p->E());                    
        }

      }

      Out_Pdg.push_back(Aux_Pdg);
      Out_PX .push_back(Aux_PX);
      Out_PY .push_back(Aux_PY);
      Out_PZ .push_back(Aux_PZ);
      Out_E  .push_back(Aux_E);                    


      NInt++;

      delete event;

    }
    else LOG("VertexGenerator", pDEBUG) << "Neutrino exiting!!!";

    if ( SecNu.size()==0 ) {
      LOG("VertexGenerator", pDEBUG) << "No more secondary neutrinos!!!";
      outflux->Fill();
      SctID.clear();
      IntID.clear();
      Tgt.clear();
      X.clear();
      Y.clear();
      Z.clear();
      T.clear();
      In_Pdg.clear();
      In_PX.clear(); In_PY.clear(); In_PZ.clear(); In_E.clear();
      Out_Pdg.clear();
      Out_PX.clear(); Out_PY.clear(); Out_PZ.clear(); Out_E.clear();
      NInt = 0;
      NNu++;
    }

  } 

  TNamed * tgeo  = new TNamed( "gOptGeometry",  gOptGeometry  );
  TNamed * tdetr = new TNamed( "gOptDetRadius", to_string(gOptDetRadius) );
  TNamed * tdeth = new TNamed( "gOptDetHeight", to_string(gOptDetHeight) );
  TNamed * tdetx = new TNamed( "gOptDetPosX",   to_string(gOptDetPos[0]) );
  TNamed * tdety = new TNamed( "gOptDetPosY",   to_string(gOptDetPos[1]) );
  TNamed * tdetz = new TNamed( "gOptDetPosZ",   to_string(gOptDetPos[2]) );
  tgeo->Write();  
  tdetr->Write();  
  tdeth->Write();  
  tdetx->Write();  
  tdety->Write();  
  tdetz->Write();  
  outflux->Write();
  outfile->Close();

  // clean-up
  delete flx_driver;
  delete trj_driver;

  return 0;
}

//**************************************************************************
//**************************************************************************
//THIS FUNCTION CALCULATE DISTANCE TO DETECTOR
//**************************************************************************
//**************************************************************************
double DistToBorder (double vx, double vy, double vz, double dx, double dy, double dz, double radius, double lw, double up) 
{

  if( TMath::Abs(dz)==1 ){
    if      ( TMath::Sqrt(vx*vx+vy*vy) > radius ) { LOG("VertexGenerator", pWARN) << "Neutrino not crossing detector!"; return 0.; }
    else if ( dz== 1 && vz>up ) { LOG("VertexGenerator", pWARN) << "Detector behind neutrino!"; return 0.; }
    else if ( dz==-1 && vz<lw ) { LOG("VertexGenerator", pWARN) << "Detector behind neutrino!"; return 0.; }
    else return (dz==1) ? TMath::Abs(up-vz) : TMath::Abs(vz-lw);
  }
  else {
    double a     = dx*dx+dy*dy;
    double b     = 2.*(vx*dx+vy*dy);
    double c     = vx*vx+vy*vy-radius*radius;
    double delta = b*b-4.*a*c;
    if ( delta<0. ) { LOG("VertexGenerator", pWARN) << "Neutrino not crossing detector!"; return 0.; }
    else {
      delta = TMath::Sqrt(delta);
      double t1 = (-b-delta)/2./a;
      double t2 = (-b+delta)/2./a;
      if ( t2<0. ) { LOG("VertexGenerator", pWARN) << "Detector behind neutrino!"; return 0.; }
      double z1 = vz+dz*t1;      
      double z2 = vz+dz*t2;
      if (z2>up) { 
        if (z1>up) { LOG("VertexGenerator", pWARN) << "Neutrino not crossing detector!"; return 0.; }
        else       { assert(dz>0.); return (up-vz)/TMath::Abs(dz); }
      }
      else if (z2<lw) {
        if (z1<lw) { LOG("VertexGenerator", pWARN) << "Neutrino not crossing detector!"; return 0.; }
        else       { assert(dz<0.); return (vz-lw)/TMath::Abs(dz); }
      }
      else return t2;
    }
  }

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
        LOG("VertexGenerator", pDEBUG) << "Reading number of events to generate";
        gOptNev = atof(argv[i]);
      }          
      if(opt.compare("--geometry")==0){   
        i++;
        LOG("VertexGenerator", pDEBUG) << "Reading geometry";
        gOptGeometry = argv[i];
      }          
      if(opt.compare("--geometry-limit")==0){   
        i++;
        LOG("VertexGenerator", pDEBUG) << "Reading geometry limit";
        gOptGeoLimit = atof(argv[i]);
      }          
      if(opt.compare("--seed")==0){   
        i++;
        LOG("VertexGenerator", pDEBUG) << "Reading seed";
        gOptRanSeed = atoi(argv[i]);
      }          
      if(opt.compare("--output")==0){   
        i++;
        LOG("VertexGenerator", pDEBUG) << "Reading output name";
        gOptOutName = argv[i];
      }          
      if(opt.compare("--probe")==0){   
        i++;
        LOG("VertexGenerator", pINFO) << "Reading neutrino flavor";
        gOptPdg = atoi(argv[i]);
      }          
      if(opt.compare("--alpha")==0){   
        i++;
        LOG("VertexGenerator", pINFO) << "Reading neutrino spectrum alpha";
        gOptAlpha = atof(argv[i]);
      }          
      if(opt.compare("--costheta")==0){   
        i++;
        LOG("VertexGenerator", pINFO) << "Reading neutrino angle";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if      (saux.size()==1) {
          gOptCthmin = gOptCthmax = atof(saux[0].c_str());
        }
        else if (saux.size()==2) {
          gOptCthmin = atof(saux[0].c_str());
          gOptCthmax = atof(saux[1].c_str());
        }
        else {
          LOG("VertexGenerator", pFATAL) << "Wrong costheta argument " << argv[i] <<"; exit";
          exit(1);
        }
      }          
      if(opt.compare("--energy")==0){   
        i++;
        LOG("VertexGenerator", pINFO) << "Reading neutrino min energy";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if      (saux.size()==1) {
          gOptEmin = gOptEmax = atof(saux[0].c_str());
        }
        else if (saux.size()==2) {
          gOptEmin = atof(saux[0].c_str());
          gOptEmax = atof(saux[1].c_str());
        }
        else {
          LOG("VertexGenerator", pFATAL) << "Wrong energy argument " << argv[i] << "; exit";
          exit(1);
        }
      }          
      if(opt.compare("--cross-sections")==0){ 
        i++;
        LOG("VertexGenerator", pINFO) << "Reading cross-section file";
        gOptInpXSecFile = argv[i];  
      }
      if(opt.compare("--detector-position")==0){
        i++; 
        LOG("VertexGenerator", pINFO) << "Reading detector pos";
        std::vector<string> saux = utils::str::Split(argv[i], ",");
        if (saux.size()==3) {
          gOptDetPos[0] = atof(saux[0].c_str());
          gOptDetPos[1] = atof(saux[1].c_str());
          gOptDetPos[2] = atof(saux[2].c_str());
        }
        else {
          LOG("VertexGenerator", pFATAL) << "Wrong detector argument " << argv[i] << "; exit";
          exit(1);
        }
      }       
      if(opt.compare("--detector-radius")==0){
        i++; 
        LOG("VertexGenerator", pINFO) << "Reading detector radius";
        gOptDetRadius = atof(argv[i]);  
      }
      if(opt.compare("--detector-height")==0){
        i++; 
        LOG("VertexGenerator", pINFO) << "Reading detector height";
        gOptDetHeight = atof(argv[i]);  
      }       
      if(opt.compare("--offset")==0){
        i++; 
        LOG("VertexGenerator", pINFO) << "Reading offset";
        gOptOffset = atof(argv[i]);  
      }       

    }
  }



  if ( gSystem->AccessPathName(gOptGeometry.c_str()) ) {
    LOG("VertexGenerator", pFATAL) << "Geometry file not found: " << gOptGeometry << "; exit";
    exit(1);
  }

  if ( gSystem->AccessPathName(gOptInpXSecFile.c_str()) ) {
    LOG("VertexGenerator", pFATAL) << "Cross section file not found: " << gOptInpXSecFile << "; exit";
    exit(1);
  }

  if ( !pdg::IsNeutrino(abs(gOptPdg)) && gOptPdg!=-1 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong neutrino flavor: " << gOptPdg << "; exit";
    exit(1);
  }

  if ( gOptNev<=0 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong number of events: " << gOptNev << "; exit";
    exit(1);    
  }

  if ( gOptGeoLimit<=0 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong geomtry limit: " << gOptGeoLimit << "; exit";
    exit(1);    
  }

  if ( gOptDetHeight<0 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong height: " << gOptDetHeight << "; exit";
    exit(1);    
  }

  if ( gOptDetRadius<0 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong radius: " << gOptDetRadius << "; exit";
    exit(1);    
  }

  if ( gOptOffset<0 ) {
    LOG("VertexGenerator", pFATAL) << "Wrong offset: " << gOptOffset << "; exit";
    exit(1);    
  }

  if (gOptEmin>gOptEmax || gOptEmin<spline_Erange.min || gOptEmax>spline_Erange.max) {
    LOG("VertexGenerator", pFATAL) << "Wrong energy range: " << gOptEmin << "," << gOptEmax << "; exit";
    exit(1);    
  }
  if (gOptEmin==gOptEmax) LOG("VertexGenerator", pNOTICE) << "Running in monoenergetic mode";

  if (gOptCthmin>gOptCthmax || gOptCthmin<-1. || gOptCthmax<-1. || gOptCthmin>1. || gOptCthmax>1. ) {
    LOG("VertexGenerator", pFATAL) << "Wrong costheta range: " << gOptCthmin << "," << gOptCthmax << "; exit";
    exit(1);    
  }
  if (gOptCthmin==gOptCthmax) LOG("VertexGenerator", pNOTICE) << "Running in monoangular mode";


  LOG("VertexGenerator", pNOTICE) << "\n\n" << utils::print::PrintFramedMesg("NuPropEarth job configuration");
  LOG("VertexGenerator", pNOTICE) 
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
     << "\n\n";


}



