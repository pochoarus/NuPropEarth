#include <cassert>
#include <cstdlib>
#include <cctype>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
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

#include "Tauola/Tauola.h"
#include "Tauola/TauolaHEPEVTParticle.h"

#include "Driver/GTRJDriver.h"
#include "Flux/IncomingFlux.h"

#include <boost/preprocessor/stringize.hpp>
extern "C" { 
  void initialize_tausic_sr_(int *muin);
  void tau_transport_sr_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
}

using namespace Tauolapp;

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::ostringstream;

using namespace genie;
using namespace genie::flux;
using namespace genie::geometry;

double d_lifetime = Tauola::tau_lifetime * 1e-3; //lifetime from mm(tauola) to m
double cSpeed = constants::kLightSpeed/(units::meter/units::second);

//functions
void GetCommandLineArgs (int argc, char ** argv);
void BuildEarth     (string geofilename);

// User-specified options:
string          gOptOutDir = "./";               
long int        gOptRanSeed;               
string          gOptInpXSecFile;           
int             gOptNev;              
int             gOptPdg;                   
double          gOptEmono = -1; //no monoenergetic
double          gOptEmin = 1e2;
double          gOptEmax = 1e10;
double          gOptCthmono = -1; //no fixed angle                   
double          gOptCthmin = -1.;
double          gOptCthmax =  1.;
double          gOptAlpha = 1; //flat spectrum in log10e
bool            gOptEnableEnergyLoss  = false;
bool            gOptEnableDecayLength = false;



const int NMAXINT = 100;

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

  string geofilename = "geometry-Earth.root";
  if ( gSystem->AccessPathName(geofilename.c_str()) ) {
    LOG("ComputeAttenuation", pDEBUG) << "Building Geometry...";
    BuildEarth(geofilename); 
  }
  geometry::ROOTGeomAnalyzer * rgeom = new geometry::ROOTGeomAnalyzer(geofilename);
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

  LOG("ComputeAttenuation", pDEBUG) << "Creating GFluxI...";
  IncomingFlux * flx_driver = new IncomingFlux(gOptPdg, gOptAlpha, gOptCthmin, gOptCthmax, gOptCthmono, gOptEmin, gOptEmax, gOptEmono);

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
  trj_driver->Configure();

  LOG("ComputeAttenuation", pDEBUG) << "Initializing TAUOLA...";
  Tauola::initialize();
  double mtau = Tauola::getTauMass(); //use this mass to avoid energy conservation warning in tauola

  LOG("ComputeAttenuation", pDEBUG) << "Initializing TAUSIC...";
  string cp_command("cp ./tausic/2004/tausic*.dat .");
  system(cp_command.c_str());

  int tauin =2;
  initialize_tausic_sr_(&tauin);

  int TAUMODEL = 1;
  int ITFLAG = 0; //set lifetime inside tausic
  double RHOR = 2.65; //rock density

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Random...";
  RandomGen * rnd = RandomGen::Instance();

  LOG("ComputeAttenuation", pINFO) << "Creating output name";
  string OutEvFile = gOptOutDir + "/NuEarthProp_" + RunOpt::Instance()->Tune()->Name() + Form("_nu%d",gOptPdg);
  if (gOptCthmono!=-1.) OutEvFile += Form("_cth%g",gOptCthmono);
  else                  OutEvFile += Form("_cth%g-%g",gOptCthmin,gOptCthmax);
  if (gOptEmono!=-1.)   OutEvFile += Form("_e%g",gOptEmono);
  else                  OutEvFile += Form("_e%g-%g",gOptEmin,gOptEmax);
  OutEvFile += Form("_s%d",gOptRanSeed);
  OutEvFile += ".root";

  LOG("ComputeAttenuation", pNOTICE) << "@@ Output file name: " << OutEvFile;

  // create file so tree is saved inside
  TFile * outfile = new TFile(OutEvFile.c_str(),"RECREATE");

  int Nu_In;
  double X4_In[4];
  double P4_In[4];
  int Nu_Out;
  double X4_Out[4];
  double P4_Out[4];
  int NInt;
  int NIntCC;
  int NIntNC;
  int SctID[NMAXINT];
  int IntID[NMAXINT];
  int Nu_pdg[NMAXINT];
  int Tgt_pdg[NMAXINT];
  double Nu_X4[NMAXINT][4];
  double Nu_P4[NMAXINT][4];
  TTree *outtree = new TTree("Events","Events");
  outtree->Branch("Nu_In",    &Nu_In,    "Nu_In/I"       );       
  outtree->Branch("P4_In",     P4_In,    "P4_In[4]/D"    );       
  outtree->Branch("X4_In",     X4_In,    "X4_In[4]/D"    );       
  outtree->Branch("Nu_Out",   &Nu_Out,   "Nu_Out/I"      );       
  outtree->Branch("P4_Out",    P4_Out,   "P4_Out[4]/D"   );       
  outtree->Branch("X4_Out",    X4_Out,   "X4_Out[4]/D"   );       
  outtree->Branch("NInt",     &NInt,     "NInt/I"        );       
  outtree->Branch("NIntCC",   &NIntCC,   "NIntCC/I"      );       
  outtree->Branch("NIntNC",   &NIntNC,   "NIntNC/I"      );       
  outtree->Branch("SctID",     SctID,    "SctID[NInt]/I" );       
  outtree->Branch("IntID",     IntID,    "IntID[NInt]/I" );       
  outtree->Branch("Nu",        Nu_pdg,   "Nu[NInt]/I"    );       
  outtree->Branch("Tgt",       Tgt_pdg,  "Tgt[NInt]/I"   );        
  outtree->Branch("X4",        Nu_X4,    "X4[NInt][4]/D" );       
  outtree->Branch("P4",        Nu_P4,    "P4[NInt][4]/D" );

  LOG("ComputeAttenuation", pDEBUG) << "Event loop...";

  int NNu = 0;
  NIntCC = 0;
  NIntNC = 0;
  NInt = 0;
  while ( NNu<gOptNev ) {

    EventRecord * event = trj_driver->GenerateEvent();

    if (NInt==0) {
      if ( NNu%100==0 ) LOG("ComputeAttenuation", pINFO) << "Event " << NNu << " out of " << gOptNev;
      Nu_In    = flx_driver->GetPdgI();
      X4_In[0] = flx_driver->GetX4I().X();  X4_In[1] = flx_driver->GetX4I().Y();  X4_In[2] = flx_driver->GetX4I().Z();  X4_In[3] = flx_driver->GetX4I().T();
      P4_In[0] = flx_driver->GetP4I().Px(); P4_In[1] = flx_driver->GetP4I().Py(); P4_In[2] = flx_driver->GetP4I().Pz(); P4_In[3] = flx_driver->GetP4I().E();
    }

    LOG("ComputeAttenuation", pDEBUG) << "Shooting Neutrino: " << NNu << "( " << NInt << " ) ----->" ;
    LOG("ComputeAttenuation", pDEBUG) << "@@Flux = " << flx_driver->GetPdgI() << " , E =" << flx_driver->GetP4I().E();
    LOG("ComputeAttenuation", pDEBUG) << "Postion   = [ " << flx_driver->GetX4I().X() << " , " << flx_driver->GetX4I().Y() << " , " << flx_driver->GetX4I().Z() << " , " << flx_driver->GetX4I().T() << " ] ";
    LOG("ComputeAttenuation", pDEBUG) << "Direction = [ " << flx_driver->GetP4I().Px()/flx_driver->GetP4I().E() << " , " << flx_driver->GetP4I().Py()/flx_driver->GetP4I().E() << " , " << flx_driver->GetP4I().Pz()/flx_driver->GetP4I().E() << " ]";

    if (event) {

      int sctid = event->Summary()->ProcInfoPtr()->ScatteringTypeId();
      int intid = event->Summary()->ProcInfoPtr()->InteractionTypeId();
      int nupdg = event->Probe()->Pdg();
      int tgt = ( event->TargetNucleus() ) ? event->TargetNucleus()->Pdg() : event->HitNucleon()->Pdg();
      const TLorentzVector & X4  = *event->Vertex();
      const TLorentzVector & P4  = *event->Probe()->P4();

      // print-out
      LOG("ComputeAttenuation", pDEBUG) << "Neutrino interaction!!!";
      LOG("ComputeAttenuation", pDEBUG) << "@@Interact   = "   << sctid << "   " << intid;
      LOG("ComputeAttenuation", pDEBUG) << "@@Probe      = "   << nupdg << " , E =" << P4.E();
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << X4.X() << " , " << X4.Y() << " , " << X4.Z() << " , " << X4.T() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << P4.Px()/P4.E() << " , " << P4.Py()/P4.E() << " , " << P4.Pz()/P4.E() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << "@@Target     = "   << tgt;
      LOG("ComputeAttenuation", pDEBUG) << *event;

      SctID[NInt]     = sctid;
      IntID[NInt]     = intid;
      Nu_pdg[NInt]    = nupdg;
      Tgt_pdg[NInt]   = tgt;
      Nu_X4[NInt][0]  = X4.X();  Nu_X4[NInt][1] = X4.Y();  Nu_X4[NInt][2] = X4.Z();  Nu_X4[NInt][3] = X4.T(); 
      Nu_P4[NInt][0]  = P4.Px(); Nu_P4[NInt][1] = P4.Py(); Nu_P4[NInt][2] = P4.Pz(); Nu_P4[NInt][3] = P4.E(); 
      
      if (intid==2) NIntCC++;
      else          NIntNC++;
      NInt++;

      int    SecNu_Pdg = 0;
      double SecNu_E = -1.; //setting energy to avoid problems when looping 
      double SecNu_Pos[4] = { 0, 0, 0, 0 };
      double SecNu_Mom[3] = { 0, 0, 0 };
      
      if (TMath::Abs(nupdg)==16 && intid==2) {

        TObjArrayIter piter(event);
        GHepParticle * p = 0;
        while ( (p=(GHepParticle *)piter.Next()) ) {
          if ( TMath::Abs(p->Pdg())!=15 ) continue; //we want the outgoing tau
          if ( p->Status()!=1           ) continue; //not final state
          if ( p->FirstMother()!=0      ) continue; 

          if (gOptEnableEnergyLoss) {
            double vxi = p->Vx()*1e-13; //from fm(genie) to cm(tausic)
            double vyi = p->Vy()*1e-13;
            double vzi = p->Vz()*1e-13;
            double ti  = p->Vt()*1e9; //from s(genie) to ns(tausic)

            const TLorentzVector & P4tau  = *p->P4();
            double dxi = P4tau.Px()/P4tau.P();
            double dyi = P4tau.Py()/P4tau.P();
            double dzi = P4tau.Pz()/P4tau.P();
            double Ei  = TMath::Sqrt( P4tau.P()*P4tau.P() + mtau*mtau );

            double depthi = 0.;
            std::vector< std::pair<double, const TGeoMaterial*> > MatLengths = geom_driver->ComputeMatLengths(X4,P4tau);
            for ( auto sitr = MatLengths.begin(); sitr != MatLengths.end(); ++sitr) {
              double length            = sitr->first * 1e2;  //from m(geom) to cm(tausic)
              const  TGeoMaterial* mat = sitr->second;
              double rho               = mat->GetDensity();
              depthi += length*rho;
            }
            depthi = depthi/RHOR;
            LOG("ComputeAttenuation", pDEBUG) << "PathLength = " << depthi << " cm r.e." ;

            double vxf,vyf,vzf,tf;
            double dxf,dyf,dzf,Ef;
            double depthf; //total distance travel after propagation

            int idec; //decay flag (0=not decay // >0=decay)
            double tauti,tautf; //dummy only used when ITFLAG!=0

            LOG("ComputeAttenuation", pDEBUG) << "Before -> X4 = " << vxi << "  " << vyi << "  " << vzi << "  " << ti;
            LOG("ComputeAttenuation", pDEBUG) << "Before -> P4 = " << dxi << "  " << dyi << "  " << dzi << "  " << Ei;

            tau_transport_sr_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&Ei,&depthi,&ti,&RHOR,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&Ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);

            LOG("ComputeAttenuation", pDEBUG) << "After -> X4 = " << vxf/100. << "  " << vyf/100. << "  " << vzf/100. << "  " << tf/1e9;
            LOG("ComputeAttenuation", pDEBUG) << "After -> P4 = " << dxf << "  " << dyf << "  " << dzf << "  " << Ef;

            if (idec>0) {
              double momf = TMath::Sqrt( Ef*Ef - mtau*mtau );
              TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
              TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( p->Pdg(), 1, dxf*momf, dyf*momf, dzf*momf, Ef, mtau, -1, -1, -1, -1 );
              Tauola_evt->addParticle(Tauola_tau);
              double pol = 1; //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped
              Tauola::decayOne(Tauola_tau,true,0.,0.,pol);
              SecNu_Pdg    = Tauola_evt->getParticle(1)->getPdgID();
              SecNu_Mom[0] = Tauola_evt->getParticle(1)->getPx(); SecNu_Mom[1] = Tauola_evt->getParticle(1)->getPy(); SecNu_Mom[2] = Tauola_evt->getParticle(1)->getPz(); 
              SecNu_E      = Tauola_evt->getParticle(1)->getE();
              SecNu_Pos[0] = vxf/100.; //from cm (tausic) to m
              SecNu_Pos[1] = vyf/100.;
              SecNu_Pos[2] = vzf/100.;
              SecNu_Pos[3] = tf/1e9; //from ns (tausic) to s
              std::cout << SecNu_E/Ef << std::endl;
              delete Tauola_evt;
              break;
            }
            else {
              LOG("ComputeAttenuation", pWARN) << "Tau did not decay!!!";
              LOG("ComputeAttenuation", pWARN) << "  Eneryi      = " << Ei;
              LOG("ComputeAttenuation", pWARN) << "  Positioni   = [ " << X4.X() << " , " << X4.Y() << " , " << X4.Z() << " , " << X4.T() << " ]";
              LOG("ComputeAttenuation", pWARN) << "  Directioni  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";
              LOG("ComputeAttenuation", pWARN) << "  PathLength = " << depthi << " cm r.e." ;
              LOG("ComputeAttenuation", pWARN) << "  Eneryf      = " << Ef;
              LOG("ComputeAttenuation", pWARN) << "  Positionf   = [ " << vxf/100. << " , " << vyf/100. << " , " << vzf/100. << " , " << tf/1e9 << " ]";
              LOG("ComputeAttenuation", pWARN) << "  Directionf  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";
            }

          }
          else {

            //decay tau
            TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
            double Etau = TMath::Sqrt( p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz() + mtau*mtau ); //use this energy to avoid energy conservation in tauola
            TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( p->Pdg(), 1, p->Px(), p->Py(), p->Pz(), Etau, mtau, -1, -1, -1, -1 );
            Tauola_evt->addParticle(Tauola_tau);
            double pol = 1; //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped
            Tauola::decayOne(Tauola_tau,true,0.,0.,pol);
            SecNu_Pdg    = Tauola_evt->getParticle(1)->getPdgID();
            SecNu_Mom[0] = Tauola_evt->getParticle(1)->getPx(); SecNu_Mom[1] = Tauola_evt->getParticle(1)->getPy(); SecNu_Mom[2] = Tauola_evt->getParticle(1)->getPz(); 
            SecNu_E      = Tauola_evt->getParticle(1)->getE();
            // position based on lifetime
            double d_r = (gOptEnableDecayLength) ? -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime : 0; //m
            SecNu_Pos[0] = p->Vx()*1e-15 + p->Px()*d_r/mtau; // initial position of tau wrt vertex changing from fm(genie) to m
            SecNu_Pos[1] = p->Vy()*1e-15 + p->Py()*d_r/mtau;
            SecNu_Pos[2] = p->Vz()*1e-15 + p->Pz()*d_r/mtau;
            SecNu_Pos[3] = p->Vt()       + Etau   *d_r/mtau/cSpeed;
            delete Tauola_evt;
            break;

          }


        }

      }
      else {
        TObjArrayIter piter(event);
        GHepParticle * p = 0;
        int mother = (intid==3) ? 0 : 4; //NC(CC): the outgoing neutrino is the daughter of the incoming neutrino (resonance W boson)
        while ( (p=(GHepParticle *)piter.Next()) ) {
          if ( p->Pdg()!=nupdg          ) continue; //we do not want flavor change
          if ( p->Status()!=1           ) continue; //not final state
          if ( p->FirstMother()!=mother ) continue; 
          SecNu_Pdg    = p->Pdg();
          SecNu_Mom[0] = p->Px(); SecNu_Mom[1] = p->Py(); SecNu_Mom[2] = p->Pz(); 
          SecNu_E      = TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()); //not using E to avoid problems with neutrinos from pythia not on shell
          //position -> passing from fm [genie] to m
          SecNu_Pos[0] = p->Vx()*1e-15; SecNu_Pos[1] = p->Vy()*1e-15; SecNu_Pos[2] = p->Vz()*1e-15; SecNu_Pos[3] = p->Vt(); 
          break;
        }
      }
      delete event;

      if ( SecNu_Pdg!=0 && SecNu_E>gOptEmin ) {
        LOG("ComputeAttenuation", pDEBUG) << "@@SecNu      = "   << SecNu_Pdg << " , E = " << SecNu_E;
        LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << SecNu_Pos[0] << " , " << SecNu_Pos[1] << " , " << SecNu_Pos[2] << " , " << SecNu_Pos[3] << " ]";
        LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << SecNu_Mom[0]/SecNu_E << " , " << SecNu_Mom[1]/SecNu_E << " , " << SecNu_Mom[2]/SecNu_E << " ]";
        flx_driver->InitNeutrino(SecNu_Mom[0],SecNu_Mom[1],SecNu_Mom[2],SecNu_E,SecNu_Pos[0],SecNu_Pos[1],SecNu_Pos[2],SecNu_Pos[3],SecNu_Pdg,X4);
        continue;
      } 

      LOG("ComputeAttenuation", pDEBUG) << " ----> Absorbed by the Earth";
      Nu_Out    = 0;
      P4_Out[0] = 0; P4_Out[1] = 0; P4_Out[2] = 0; P4_Out[3] = 0;
      X4_Out[0] = 0; X4_Out[1] = 0; X4_Out[2] = 0; X4_Out[3] = 0;

    }
    else {
      LOG("ComputeAttenuation", pDEBUG) << " ----> Goodbye Earth!!!";
      P4_Out[0] = flx_driver->GetP4I().Px(); P4_Out[1] = flx_driver->GetP4I().Py(); P4_Out[2] = flx_driver->GetP4I().Pz(); P4_Out[3] = flx_driver->GetP4I().E();
      X4_Out[0] = flx_driver->GetX4I().X(); X4_Out[1] = flx_driver->GetX4I().Y(); X4_Out[2] = flx_driver->GetX4I().Z(); X4_Out[3] = flx_driver->GetX4I().T();
      Nu_Out = (P4_Out[3]==P4_In[3] && X4_Out[2]==X4_In[2]) ? 1 : flx_driver->GetPdgI();     
    }

    outtree->Fill();

    NInt = 0;
    NIntCC = 0;
    NIntNC = 0;
    NNu++;

  } 

  outtree->Write();
  outfile->Close();

  // clean-up
  delete geom_driver;
  delete flx_driver;
  delete trj_driver;

  system("rm geometry-Earth.root");
  system("rm tausic*.dat");

  return 0;
}

//**************************************************************************
//**************************************************************************
//THIS FUNCTION BUILD THE GEOMETRY OF THE EARTH
//**************************************************************************
//**************************************************************************
void BuildEarth(string geofilename){

  double fREarth_km = fREarth_m/1e3;

  LOG("ComputeAttenuation", pINFO) << "fREarth_km = "<< fREarth_km;

  struct Layer{ double r1, r2, rho; TString Composition; };
  vector<Layer> VecLayers;

  const int NLAYERS = 9;
  TString fComp[NLAYERS] = { "Core", "Core", "Mantle", "Mantle", "Mantle", "Mantle", "Mantle", "Mantle", "Rock"       };
  double fR[NLAYERS]     = { 1221.5,  3480.,    5701.,    5771.,    5971.,    6151.,   6346.6,    6356.,  fREarth_km }; //km

  double fEarthCoeff[NLAYERS][4] = {
    {13.0885,0.,-8.8381,0.},
    {12.5815,-1.2638,-3.6426,-5.5281} ,
    {7.9565,-6.4761,5.5283,-3.0807},
    {5.3197,-1.4836,0.,0.},
    {11.2494,-8.0298,0.,0.},
    {7.1089,-3.8045,0.,0.},
    {2.691,0.6924,0.,0.},
    {2.9,0.,0.,0.},
    {2.65,0.,0.,0.}
  };
  
  
  double rmean = 0.;
  double step  = 100.;
  double r1    = 0;
  double r2    = step;  
  bool border  = false;  
  int iLayer   = 0;
  int iStep    = 0;
    
  while(1){       

    if( r2>fR[iLayer] ) {
      r2     = fR[iLayer]; 
      border = true;
    }
   
    Layer layer;
    layer.r1 = r1;
    layer.r2 = r2;
    rmean    = ( r2 + r1 ) / 2.;
    layer.rho         = fEarthCoeff[iLayer][0] + fEarthCoeff[iLayer][1]*rmean/fREarth_km + fEarthCoeff[iLayer][2]*pow(rmean/fREarth_km,2) + fEarthCoeff[iLayer][3]*pow(rmean/fREarth_km,3);
    layer.Composition = fComp[iLayer];
    
    VecLayers.push_back(layer);
    
    if( border ) {
      r1     = r2;
      border = false;
      iLayer++;
    }
    else r1 = (iStep+1)*step;

    r2 = (iStep+2)*step;
    iStep++;
    
    if( r1>=fR[NLAYERS-3] ) break;
    
  }

  for(int iiLayer=NLAYERS-2; iiLayer<NLAYERS; iiLayer++){ // constant density
    Layer layer;
    layer.r1          = fR[iiLayer-1];
    layer.r2          = fR[iiLayer];
    layer.rho         = fEarthCoeff[iiLayer][0] + fEarthCoeff[iiLayer][1]*layer.r2/fREarth_km + fEarthCoeff[iiLayer][2]*pow(layer.r2/fREarth_km,2) + fEarthCoeff[iiLayer][3]*pow(layer.r2/fREarth_km,3);
    layer.Composition = fComp[iiLayer];    
    VecLayers.push_back(layer);
  }

  map<int,double> RockComp, MantleComp, CoreComp;
  RockComp    .insert(map<int, double>::value_type(1000080160,0.463)     );
  RockComp    .insert(map<int, double>::value_type(1000140280,0.282)     );
  RockComp    .insert(map<int, double>::value_type(1000130270,0.0823)    );
  RockComp    .insert(map<int, double>::value_type(1000260560,0.0563)    );
  RockComp    .insert(map<int, double>::value_type(1000200400,0.0415)    );
  RockComp    .insert(map<int, double>::value_type(1000110230,0.0236)    );
  RockComp    .insert(map<int, double>::value_type(1000120240,0.0233)    );
  RockComp    .insert(map<int, double>::value_type(1000190390,0.0209)    );
  RockComp    .insert(map<int, double>::value_type(1000220480,0.0057)    );
  RockComp    .insert(map<int, double>::value_type(1000010010,0.0014)    );
  MantleComp  .insert(map<int, double>::value_type(1000080160,0.4522)    );
  MantleComp  .insert(map<int, double>::value_type(1000120240,0.2283)    );
  MantleComp  .insert(map<int, double>::value_type(1000140280,0.2149)    );
  MantleComp  .insert(map<int, double>::value_type(1000260560,0.0597)    );
  MantleComp  .insert(map<int, double>::value_type(1000130270,0.0225)    );
  MantleComp  .insert(map<int, double>::value_type(1000200400,0.0224)    );
  CoreComp    .insert(map<int, double>::value_type(1000260560,0.9)       );
  CoreComp    .insert(map<int, double>::value_type(1000280580,0.1)       );

  // Define media and volumes
  TGeoManager * GeoManager =  new TGeoManager("VolGenGeo", "generation volume geometry");

  TGeoTranslation * Trans = new TGeoTranslation(0.,0.,0.);
  
  // vacuum
  TGeoMaterial * MatVacuum = new TGeoMaterial("MatVacuum");
  MatVacuum->SetA(0.);
  MatVacuum->SetZ(0.);
  MatVacuum->SetDensity(0.);
  TGeoMedium * Vacuum = new TGeoMedium("Vacuum", 0, MatVacuum);
  TGeoVolume *TopVolume = GeoManager->MakeTube("TopVolume",Vacuum,0.,fR[NLAYERS-1]*1E3, 2*fR[NLAYERS-1]*1.E3);
  GeoManager->SetTopVolume(TopVolume);

  TGeoMixture * LayerMix    [VecLayers.size()]; 
  TGeoMedium  * LayerMedium [VecLayers.size()];
  TGeoVolume  * Layer       [VecLayers.size()];
  for(unsigned iiLayer=0; iiLayer<VecLayers.size(); iiLayer++){
    
    LOG("ComputeAttenuation", pDEBUG) << "Layer=" << iiLayer << ", Composition=" << VecLayers[iiLayer].Composition << ", r1=" << VecLayers[iiLayer].r1 << ", r2=" << VecLayers[iiLayer].r2 << ", rho=" << VecLayers[iiLayer].rho;

    TString name = VecLayers[iiLayer].Composition;
    name += iiLayer;
    
    map<int,double>::iterator iter;
    if ( VecLayers[iiLayer].Composition=="Rock" )     {
      LayerMix[iiLayer] = new TGeoMixture( name, RockComp.size(), VecLayers[iiLayer].rho ); 
      iter = RockComp.begin();
      for( ; iter != RockComp.end(); ++iter )     LayerMix[iiLayer]->AddElement( pdg::IonPdgCodeToA(iter->first), pdg::IonPdgCodeToZ(iter->first), iter->second );                  
    }
    else if( VecLayers[iiLayer].Composition=="Mantle" )    {
      LayerMix[iiLayer] = new TGeoMixture( name, MantleComp.size(), VecLayers[iiLayer].rho ); 
      iter = MantleComp.begin();
      for( ; iter != MantleComp.end(); ++iter )   LayerMix[iiLayer]->AddElement( pdg::IonPdgCodeToA(iter->first), pdg::IonPdgCodeToZ(iter->first), iter->second );                  
    }
    else if( VecLayers[iiLayer].Composition=="Core" )      {
      LayerMix[iiLayer] = new TGeoMixture( name, CoreComp.size(), VecLayers[iiLayer].rho ); 
      iter = CoreComp.begin();
      for( ; iter != CoreComp.end(); ++iter )     LayerMix[iiLayer]->AddElement( pdg::IonPdgCodeToA(iter->first), pdg::IonPdgCodeToZ(iter->first), iter->second );                  
    }  
    
    LayerMedium[iiLayer] = new TGeoMedium( name, iiLayer+1, LayerMix[iiLayer] );   
    Layer[iiLayer]       = GeoManager->MakeSphere( name, LayerMedium[iiLayer], VecLayers[iiLayer].r1*1.E3, VecLayers[iiLayer].r2*1.E3 ); 
    TopVolume->AddNode( Layer[iiLayer], iiLayer+1, Trans );
    
  }

  GeoManager->CloseGeometry();
  GeoManager->Export(geofilename.c_str());

  Vacuum    = 0; delete Vacuum;
  MatVacuum = 0; delete MatVacuum;
  for(unsigned iiLayer=0; iiLayer<VecLayers.size(); iiLayer++){
    LayerMedium[iiLayer] = 0; delete LayerMedium[iiLayer];
    LayerMix[iiLayer]    = 0; delete LayerMix[iiLayer];        
  }  
  GeoManager = 0; delete GeoManager;

  return;

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

      if(opt.compare("-n")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading number of events to generate";
        gOptNev = atof(argv[i]);
      }          
      if(opt.compare("-s")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading seed";
        gOptRanSeed = atoi(argv[i]);
      }          
      if(opt.compare("-o")==0){   
        i++;
        LOG("ComputeAttenuation", pDEBUG) << "Reading output dir";
        gOptOutDir = argv[i];
      }          
      if(opt.compare("-p")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino flavor";
        gOptPdg = atoi(argv[i]);
      }          
      if(opt.compare("-a")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino spectrum alpha";
        gOptAlpha = atof(argv[i]);
      }          
      if(opt.compare("-t")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino mono angle";
        gOptCthmono = atof(argv[i]);
      }          
      if(opt.compare("-tmin")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino min angle";
        gOptCthmin = atof(argv[i]);
      }          
      if(opt.compare("-tmax")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino max angle";
        gOptCthmax = atof(argv[i]);
      }          
      if(opt.compare("-e")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino mono energy";
        gOptEmono = atof(argv[i]);
      }          
      if(opt.compare("-emin")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino min energy";
        gOptEmin = atof(argv[i]);
      }          
      if(opt.compare("-emax")==0){   
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading neutrino max energy";
        gOptEmax = atof(argv[i]);
      }          
      if(opt.compare("--cross-sections")==0){ 
        i++;
        LOG("ComputeAttenuation", pINFO) << "Reading cross-section file";
        gOptInpXSecFile = argv[i];  
      } 
      if(opt.compare("-enable-eloss")==0){ 
        LOG("ComputeAttenuation", pINFO) << "Enable energy loss for tau propagation";
        gOptEnableEnergyLoss = true;  
      } 
      if(opt.compare("-enable-decaylength")==0){ 
        LOG("ComputeAttenuation", pINFO) << "Enable decay length for tau propagation";
        gOptEnableDecayLength = true;  
      } 

    }
  }

  ostringstream expinfo;
  if(gOptNev > 0)            { expinfo << gOptNev            << " events";   } 

  LOG("ComputeAttenuation", pNOTICE) << "\n\n" << utils::print::PrintFramedMesg("NuPropEarth job configuration");

  LOG("ComputeAttenuation", pNOTICE) 
     << "\n"
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Using cross-section file: " << gOptInpXSecFile
     << "\n @@ Using output dir: " << gOptOutDir
     << "\n @@ Exposure" 
     << "\n\t" << expinfo.str()
     << "\n @@ Kinematics" 
     << "\n\t Pdg          = " << gOptPdg
     << "\n\t CosTheta     = " << gOptCthmono
     << "\n\t CosTheta_min = " << gOptCthmin
     << "\n\t CosTheta_max = " << gOptCthmax
     << "\n\t E            = " << gOptEmono << " GeV"
     << "\n\t Emin         = " << gOptEmin << " GeV"
     << "\n\t Emax         = " << gOptEmax << " GeV"
     << "\n\t Alpha        = " << gOptAlpha
     << "\n\t EnergyLoss   = " << gOptEnableEnergyLoss
     << "\n\t DecayLength  = " << gOptEnableDecayLength
     << "\n\n";


}



