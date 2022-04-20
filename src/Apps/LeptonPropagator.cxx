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
#include "Propagation/MuonPropagation.h"
#include "Propagation/HadronPropagation.h"


using std::string;
using std::pair;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

//functions
void GetCommandLineArgs (int argc, char ** argv);
bool IsInDetector       (double vx, double vy, double vz, double * detpos, double radius, double lw, double up);

// User-specified options:
string          gOptInName     = "./in.root";               
string          gOptOutName    = "./out.root";               
int             gOptRanSeed    = 0;               
string          gOptTauProp    = ""; //Path to proposal tables
string          gOptMuonProp   = ""; //Path to proposal tables
bool            gOptHadronProp = true;
double          gOptPropStep   = 10.; //step in propagation [m]
double          gOptMaxProp    = 1e5; //maximum distance propagation [m]
double          gOptMinEnergy  = 10.; //minimum energy we propagate

//**************************************************************************
//**************************************************************************
//THIS IS THE MAIN FUNCTION
//**************************************************************************
//**************************************************************************
int main(int argc, char** argv)
{

  // Parse command line arguments
  LOG("LeptonPropagator", pDEBUG) << "Reading options...";
  GetCommandLineArgs(argc,argv);

  TFile * infile = new TFile(gOptInName.c_str());

  string          gInGeometry  = infile->Get("gOptGeometry")->GetTitle();               
  double          gInDetRadius = atof(infile->Get("gOptDetRadius")->GetTitle());
  double          gInDetHeight = atof(infile->Get("gOptDetHeight")->GetTitle());
  double          gInDetPos[3] = { atof(infile->Get("gOptDetPosX")->GetTitle()), 
                                   atof(infile->Get("gOptDetPosY")->GetTitle()), 
                                   atof(infile->Get("gOptDetPosZ")->GetTitle()) };

  LOG("LeptonPropagator", pNOTICE) << "gInGeometry  : " << gInGeometry;
  LOG("LeptonPropagator", pNOTICE) << "gInDetRadius : " << gInDetRadius;
  LOG("LeptonPropagator", pNOTICE) << "gInDetHeight : " << gInDetHeight;
  LOG("LeptonPropagator", pNOTICE) << "gInDetPos    : " << gInDetPos[0] << " , " << gInDetPos[1] << " , " << gInDetPos[2];

  geometry::ROOTGeomAnalyzer * rgeom = new geometry::ROOTGeomAnalyzer(gInGeometry);
  rgeom -> SetLengthUnits  ( genie::utils::units::UnitFromString("m")     );
  rgeom -> SetDensityUnits ( genie::utils::units::UnitFromString("g_cm3") );
  rgeom -> SetTopVolName   ("");

  TGeoVolume * topvol = rgeom->GetGeometry()->GetTopVolume();
   if(!topvol) {
    LOG("LeptonPropagator", pFATAL) << " ** Null top ROOT geometry volume!";
    exit(1);
  }
  LOG("LeptonPropagator", pNOTICE) << topvol->GetName();


  LOG("ComputeAttenuation", pNOTICE) << "Initializing Hadron Propagation...";
  HadronPropagation * hadronprop = new HadronPropagation(rgeom);

  LOG("ComputeAttenuation", pNOTICE) << "Initializing Tau Propagation...";
  TauPropagation * tauprop = new TauPropagation(gOptTauProp,gOptRanSeed,rgeom,{"CoreA","CoreB","MantleA","MantleB","MantleC","MantleD","MatVacuum"});

  LOG("ComputeAttenuation", pNOTICE) << "Initializing Muon Propagation...";
  MuonPropagation * muonprop = new MuonPropagation(gOptMuonProp,gOptRanSeed,rgeom,{"CoreA","CoreB","MantleA","MantleB","MantleC","MantleD","MatVacuum"});

  TDatabasePDG * PdgDB  = TDatabasePDG::Instance();

  double Wgt;
  vector<double> *X = 0;
  vector<double> *Y = 0;
  vector<double> *Z = 0;
  vector<double> *In_Pdg = 0;
  vector<double> *In_PX  = 0;
  vector<double> *In_PY  = 0;
  vector<double> *In_PZ  = 0;
  vector<double> *In_E   = 0;
  vector<vector<double>> *Out_Pdg = 0;
  vector<vector<double>> *Out_PX  = 0;
  vector<vector<double>> *Out_PY  = 0;
  vector<vector<double>> *Out_PZ  = 0;
  vector<vector<double>> *Out_E   = 0;

  TTree *outflux = (TTree*)infile->Get("OutFlux");
  outflux->SetBranchAddress("Wgt",     &Wgt     );
  outflux->SetBranchAddress("X",       &X       );
  outflux->SetBranchAddress("Y",       &Y       );
  outflux->SetBranchAddress("Z",       &Z       );
  outflux->SetBranchAddress("In_Pdg",  &In_Pdg  );
  outflux->SetBranchAddress("In_PX",   &In_PX   );
  outflux->SetBranchAddress("In_PY",   &In_PY   );
  outflux->SetBranchAddress("In_PZ",   &In_PZ   );
  outflux->SetBranchAddress("In_E",    &In_E    );
  outflux->SetBranchAddress("Out_Pdg", &Out_Pdg );
  outflux->SetBranchAddress("Out_PX",  &Out_PX  );
  outflux->SetBranchAddress("Out_PY",  &Out_PY  );
  outflux->SetBranchAddress("Out_PZ",  &Out_PZ  );
  outflux->SetBranchAddress("Out_E",   &Out_E   );

  TFile * outfile = new TFile(gOptOutName.c_str(),"RECREATE");
  outfile->cd();
  
  double PrWgt;
  int PrPdg;
  double PrPX,PrPY,PrPZ,PrE;
  vector<int>    Det_In_Pdg;
  vector<double> Det_In_X,Det_In_Y,Det_In_Z;
  vector<double> Det_In_PX,Det_In_PY,Det_In_PZ,Det_In_E;
  vector<vector<double>> Det_Out_Pdg;
  vector<vector<double>> Det_Out_X,Det_Out_Y,Det_Out_Z;
  vector<vector<double>> Det_Out_PX,Det_Out_PY,Det_Out_PZ,Det_Out_E;

  TTree * propflux = new TTree("PropFlux","PropFlux");
  propflux->Branch("PrWgt",   &PrWgt, "PrWgt/D" );
  propflux->Branch("PrPdg",   &PrPdg, "PrPdg/I" );
  propflux->Branch("PrPX",    &PrPX,  "PrPX/D"  );
  propflux->Branch("PrPY",    &PrPY,  "PrPY/D"  );
  propflux->Branch("PrPZ",    &PrPZ,  "PrPZ/D"  );
  propflux->Branch("PrE",     &PrE,   "PrE/D"   );
  propflux->Branch("In_X",    &Det_In_X     );
  propflux->Branch("In_Y",    &Det_In_Y     );
  propflux->Branch("In_Z",    &Det_In_Z     );
  propflux->Branch("In_Pdg",  &Det_In_Pdg   );
  propflux->Branch("In_PX",   &Det_In_PX    );
  propflux->Branch("In_PY",   &Det_In_PY    );
  propflux->Branch("In_PZ",   &Det_In_PZ    );
  propflux->Branch("In_E",    &Det_In_E     );
  propflux->Branch("Out_X",   &Det_Out_X    );
  propflux->Branch("Out_Y",   &Det_Out_Y    );
  propflux->Branch("Out_Z",   &Det_Out_Z    );
  propflux->Branch("Out_Pdg", &Det_Out_Pdg  );
  propflux->Branch("Out_PX",  &Det_Out_PX   );
  propflux->Branch("Out_PY",  &Det_Out_PY   );
  propflux->Branch("Out_PZ",  &Det_Out_PZ   );
  propflux->Branch("Out_E",   &Det_Out_E    );

  int Nevents = outflux->GetEntries();

  for ( int e=0; e<Nevents; e++ ) {
    
    if ( e%100==0 ) LOG("LeptonPropagator", pNOTICE) << "Event " << e << " out of " << Nevents;

    outflux->GetEntry(e);

    PrWgt = Wgt;
    PrPdg = In_Pdg->at(0);
    PrPX  = In_PX->at(0);
    PrPY  = In_PY->at(0);
    PrPZ  = In_PZ->at(0);
    PrE   = In_E->at(0);

    int NInt = X->size();

    for ( int i=0; i<NInt; i++) {

      int Nout = Out_Pdg->at(i).size();

      LOG("LeptonPropagator", pDEBUG) << "Event " << e << " (" << i << "): PrPdg = " << PrPdg << " -------------->" ;
      LOG("LeptonPropagator", pDEBUG) << "@@ In_Pdg " << In_Pdg->at(i);
      LOG("LeptonPropagator", pDEBUG) << "@@ In_E   " << In_E->at(i);
      LOG("LeptonPropagator", pDEBUG) << "@@ Nout " << Nout;

      //first lets propagate hadrons to leptons
      if (gOptHadronProp) {
        vector<GHepParticle> Sec;
        for (int o=0; o<Nout; o++) {
          if ( Out_E->at(i)[o]<gOptMinEnergy ) continue;
          GHepParticle * p = new GHepParticle(Out_Pdg->at(i)[o],kIStUndefined,-1,-1,-1,-1,Out_PX->at(i)[o],Out_PY->at(i)[o],Out_PZ->at(i)[o],Out_E->at(i)[o],X->at(i),Y->at(i),Z->at(i),0.);
          TParticlePDG * partinfo = PdgDB->GetParticle(p->Pdg());
          string pclass = partinfo->ParticleClass();
          if( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) { 
            std::vector<GHepParticle> CharmProd = hadronprop->Propagate(p,pclass);
            for ( auto & cpr : CharmProd ) {
              if      ( pdg::IsTau     (TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr);
              else if ( pdg::IsMuon    (TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr);
              else if ( pdg::IsNeutrino(TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr);
            }
          }
          else if( pclass=="B-Meson" || pclass=="B-Baryon" ) {           
            std::vector<GHepParticle> BottomProd = hadronprop->Propagate(p,pclass);
            for ( auto & bpr : BottomProd ) {
              if      ( pdg::IsTau     (TMath::Abs(bpr.Pdg())) ) Sec.push_back(bpr);
              else if ( pdg::IsMuon    (TMath::Abs(bpr.Pdg())) ) Sec.push_back(bpr);
              else if ( pdg::IsNeutrino(TMath::Abs(bpr.Pdg())) ) Sec.push_back(bpr);
              else {
                partinfo = PdgDB->GetParticle(bpr.Pdg());
                pclass = partinfo->ParticleClass();
                if( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) { 
                  std::vector<GHepParticle> CharmProd = hadronprop->Propagate(&bpr,pclass);
                  for ( auto & cpr : CharmProd ) {
                    if      ( pdg::IsTau     (TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr); 
                    else if ( pdg::IsMuon    (TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr);
                    else if ( pdg::IsNeutrino(TMath::Abs(cpr.Pdg())) ) Sec.push_back(cpr);
                  }
                }
              }
            }
          } 
          for ( auto & s : Sec ) {
            Out_Pdg->at(i).push_back(s.Pdg());
            Out_PX->at(i).push_back(s.Px());
            Out_PY->at(i).push_back(s.Py());
            Out_PZ->at(i).push_back(s.Pz());
            Out_E->at(i).push_back(s.E());
          }
          delete p;
        }
      }
      
      //now check if any neutrino is inside detector
      if ( IsInDetector(X->at(i),Y->at(i),Z->at(i),gInDetPos,gInDetRadius,-gInDetHeight/2.,gInDetHeight/2.) ) {

        LOG("LeptonPropagator", pDEBUG) << "In detector!";
        
        Det_In_Pdg.push_back(In_Pdg->at(i));
        Det_In_E.push_back(In_E->at(i));
        Det_In_PX.push_back(In_PX->at(i)); Det_In_PY.push_back(In_PY->at(i)); Det_In_PZ.push_back(In_PZ->at(i));
        Det_In_X.push_back(X->at(i)); Det_In_Y.push_back(Y->at(i)); Det_In_Z.push_back(Z->at(i));
        
        //we only store long-live leptons and neutrinos (the rest is shower like)
        vector<double> Aux_Pdg;
        vector<double> Aux_E;
        vector<double> Aux_PX,Aux_PY,Aux_PZ;
        vector<double> Aux_X,Aux_Y,Aux_Z;
        for (int o=0; o<Nout; o++) {
          if ( Out_E->at(i)[o]<gOptMinEnergy ) continue;
          int apdg = TMath::Abs(Out_Pdg->at(i)[o]); 
          if ( pdg::IsTau( apdg ) || pdg::IsMuon( apdg ) || pdg::IsNeutrino( apdg ) ) {
            LOG("LeptonPropagator", pDEBUG) << "  out: " << Out_Pdg->at(i)[o];
            Aux_Pdg.push_back(Out_Pdg->at(i)[o]);
            Aux_E.push_back(Out_E->at(i)[o]);
            Aux_PX.push_back(Out_PX->at(i)[o]); Aux_PY.push_back(Out_PY->at(i)[o]); Aux_PZ.push_back(Out_PZ->at(i)[o]);
            Aux_X.push_back(X->at(i)); Aux_Y.push_back(Y->at(i)); Aux_Z.push_back(Z->at(i));            
          }
        }
        Det_Out_Pdg.push_back(Aux_Pdg);
        Det_Out_E.push_back(Aux_E);
        Det_Out_PX.push_back(Aux_PX); Det_Out_PY.push_back(Aux_PY); Det_Out_PZ.push_back(Aux_PZ);
        Det_Out_X.push_back(Aux_X); Det_Out_Y.push_back(Aux_Y); Det_Out_Z.push_back(Aux_Z);

      }

      //finally check if any neutrino outside detector can reach detector
      else {

        vector<GHepParticle> PropPart;
        vector<GHepParticle> SecMuon;

        for (int o=0; o<Nout; o++) {

          if ( Out_E->at(i)[o]<gOptMinEnergy ) continue;

          GHepParticle * p = new GHepParticle(Out_Pdg->at(i)[o],kIStUndefined,-1,-1,-1,-1,Out_PX->at(i)[o],Out_PY->at(i)[o],Out_PZ->at(i)[o],Out_E->at(i)[o],X->at(i),Y->at(i),Z->at(i),0.);

          LOG("LeptonPropagator", pDEBUG) << "Propagating: " << p->Pdg() << "  " << p->E();
          
          if ( pdg::IsTau( TMath::Abs(p->Pdg()) ) ) {
            double length = 0;
            while ( length<gOptMaxProp ) {
              tauprop->Propagate(p,gOptPropStep,gOptMinEnergy);
              if ( IsInDetector(p->Vx(),p->Vy(),p->Vz(),gInDetPos,gInDetRadius,-gInDetHeight/2.,gInDetHeight/2.) ) {
                LOG("LeptonPropagator", pDEBUG) << "Tau reach detector: " << p->Pdg() << "  " << p->E();
                PropPart.push_back(*p);
                break; 
              }              
              if ( p->E()<gOptMinEnergy ) {
                LOG("LeptonPropagator", pDEBUG) << "Too low energy decay: " << p->E();
                break;                
              }
              if ( !pdg::IsTau( TMath::Abs(p->Pdg()) ) ) {
                LOG("LeptonPropagator", pDEBUG) << "Tau decay: " << p->Pdg();
                if ( pdg::IsMuon( TMath::Abs(p->Pdg()) ) ) SecMuon.push_back(*p);
                break;
              }
              length += gOptPropStep;
            }
          }
          else if ( pdg::IsMuon( TMath::Abs(p->Pdg()) ) ) { 
            double length = 0;
            while ( length<gOptMaxProp ) {
              muonprop->Propagate(p,gOptPropStep,gOptMinEnergy);
              if ( IsInDetector(p->Vx(),p->Vy(),p->Vz(),gInDetPos,gInDetRadius,-gInDetHeight/2.,gInDetHeight/2.) ) {
                LOG("LeptonPropagator", pDEBUG) << "Muon reach detector: " << p->Pdg() << "  " << p->E();
                PropPart.push_back(*p);
                break; 
              }
              if ( p->E()<gOptMinEnergy ) {
                LOG("LeptonPropagator", pDEBUG) << "Too low energy decay: " << p->E();
                break;                
              }
              if ( !pdg::IsMuon( TMath::Abs(p->Pdg()) ) ) {
                LOG("LeptonPropagator", pDEBUG) << "Muon decay: " << p->Pdg();
                break;
              }
              length += gOptPropStep;
            }
          }          
          delete p;
        }

        for ( auto & p : SecMuon ) {

          LOG("LeptonPropagator", pNOTICE) << "Propagating secondaries: " << p.Pdg() << " " << p.E();
          
          double length = 0;
          while ( length<gOptMaxProp ) {
            muonprop->Propagate(&p,gOptPropStep,gOptMinEnergy);
            if ( IsInDetector(p.Vx(),p.Vy(),p.Vz(),gInDetPos,gInDetRadius,-gInDetHeight/2.,gInDetHeight/2.) ) {
              LOG("LeptonPropagator", pDEBUG) << "Tau-Muon reach detector: " << p.Pdg() << " " << p.E();
              PropPart.push_back(p);
              break; 
            }
              if ( p.E()<gOptMinEnergy ) {
                LOG("LeptonPropagator", pDEBUG) << "Too low energy decay: " << p.E();
                break;                
              }
            if ( !pdg::IsMuon( TMath::Abs(p.Pdg()) ) ) {
              LOG("LeptonPropagator", pDEBUG) << "Tau-Muon decay: " << p.Pdg();
              break;
            }
            length += gOptPropStep;
          }
        }

        if ( PropPart.size()==0 ) continue;
        
        Det_In_Pdg.push_back(In_Pdg->at(i));
        Det_In_E.push_back(In_E->at(i));
        Det_In_PX.push_back(In_PX->at(i)); Det_In_PY.push_back(In_PY->at(i)); Det_In_PZ.push_back(In_PZ->at(i));
        Det_In_X.push_back(X->at(i)); Det_In_Y.push_back(Y->at(i)); Det_In_Z.push_back(Z->at(i));

        vector<double> Aux_Pdg;
        vector<double> Aux_E;
        vector<double> Aux_PX,Aux_PY,Aux_PZ;
        vector<double> Aux_X,Aux_Y,Aux_Z;
        for ( auto & p : PropPart ) {
          LOG("LeptonPropagator", pDEBUG) << "  out: " << p.Pdg();
          Aux_Pdg.push_back(p.Pdg());
          Aux_E.push_back(p.E());
          Aux_PX.push_back(p.Px()); Aux_PY.push_back(p.Py()); Aux_PZ.push_back(p.Pz());
          Aux_X.push_back(p.Vx()); Aux_Y.push_back(p.Vy()); Aux_Z.push_back(p.Vz());            
        }
        Det_Out_Pdg.push_back(Aux_Pdg);
        Det_Out_E.push_back(Aux_E);
        Det_Out_PX.push_back(Aux_PX); Det_Out_PY.push_back(Aux_PY); Det_Out_PZ.push_back(Aux_PZ);
        Det_Out_X.push_back(Aux_X); Det_Out_Y.push_back(Aux_Y); Det_Out_Z.push_back(Aux_Z);

      }

    }
    
    if ( Det_In_Pdg.size()>0 ) propflux->Fill();

    Det_In_Pdg.clear();
    Det_In_E.clear();
    Det_In_PX.clear(); Det_In_PY.clear(); Det_In_PZ.clear();
    Det_In_X.clear(); Det_In_Y.clear(); Det_In_Z.clear();
    Det_Out_Pdg.clear();
    Det_Out_E.clear();
    Det_Out_PX.clear(); Det_Out_PY.clear(); Det_Out_PZ.clear();
    Det_Out_X.clear(); Det_Out_Y.clear(); Det_Out_Z.clear();


  } 

  LOG("LeptonPropagator", pDEBUG) << "Writing output file ...";
  propflux->Write();
  outfile->Close();

  return 0;
}


//**************************************************************************
//**************************************************************************
//THIS FUNCTION CALCULATE DISTANCE TO DETECTOR
//**************************************************************************
//**************************************************************************
bool IsInDetector(double vx, double vy, double vz, double * detpos, double radius, double lw, double up) 
{
  vx -= detpos[0];
  vy -= detpos[1];
  vz -= detpos[2];
  return TMath::Sqrt(vx*vx+vy*vy)<radius && vz>lw && vz<up;
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

      if(opt.compare("--input")==0){   
        i++;
        LOG("LeptonPropagator", pDEBUG) << "Reading input file";
        gOptInName = argv[i];
      }          
      if(opt.compare("--seed")==0){   
        i++;
        LOG("LeptonPropagator", pDEBUG) << "Reading seed";
        gOptRanSeed = atoi(argv[i]);
      }          
      if(opt.compare("--output")==0){   
        i++;
        LOG("LeptonPropagator", pDEBUG) << "Reading output name";
        gOptOutName = argv[i];
      }          
      if(opt.compare("--tau-propagation")==0){ 
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading tau propagation tables";
        gOptTauProp = argv[i];  
      }      
      if(opt.compare("--muon-propagation")==0){ 
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading muon propagation tables";
        gOptMuonProp = argv[i];  
      }      
      if(opt.compare("--no-hadron-propagation")==0){ 
        LOG("ComputeAttenuation", pINFO) << "Disable hadron propagation";
        gOptHadronProp = false;  
      }      

    }
  }


  if ( gSystem->AccessPathName(gOptInName.c_str()) ) {
    LOG("LeptonPropagator", pFATAL) << "Input file not found: " << gOptInName << "; exit";
    exit(1);
  }

  if ( gSystem->AccessPathName(gOptMuonProp.c_str()) ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong muon propagation: " << gOptMuonProp << "; exit";
    exit(1);          
  }

  if ( gSystem->AccessPathName(gOptTauProp.c_str()) ) {
    LOG("ComputeAttenuation", pFATAL) << "Wrong tau propagation: " << gOptTauProp << "; exit";
    exit(1);          
  }

  LOG("LeptonPropagator", pNOTICE) << "\n\n" << utils::print::PrintFramedMesg("NuPropEarth job configuration");
  LOG("LeptonPropagator", pNOTICE) 
     << "\n"
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using output file name: " << gOptOutName
     << "\n @@ Input" 
     << "\n\t" << gOptInName
     << "\n @@ Propagation" 
     << "\n\t TauProp      = " << gOptTauProp
     << "\n\t MuonProp     = " << gOptMuonProp
     << "\n\t HadronProp   = " << gOptHadronProp
     << "\n\n";


}



