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

//TAUSIC variables
int TAUIN = 2;
int TAUMODEL = 1;
int ITFLAG = 0; //set lifetime inside tausic
double RHOR = 2.65; //rock density

//TAUOLA varible
double d_lifetime = Tauola::tau_lifetime * 1e-1; //lifetime from mm(tauola) to cm
double mtau = Tauola::getTauMass();
double cSpeed = constants::kLightSpeed/(units::centimeter/units::nanosecond);

//functions
void GetCommandLineArgs (int argc, char ** argv);
void BuildEarth     (string geofilename);
void FillNeutrino       (GHepParticle * nu, int &pdg, double &px, double &py,double &pz, double &e, double &vx, double &vy, double &vz, double &t) {
  pdg = nu->Pdg();
  px  = nu->Px(); py = nu->Py(); pz = nu->Pz(); e = nu->E();
  vx  = nu->Vx(); vy = nu->Vy(); vz = nu->Vz(); t = nu->Vt();
}
vector<GHepParticle> DecayTau (GHepParticle * tau);


// User-specified options:
string          gOptOutDir = "./";               
long int        gOptRanSeed;               
string          gOptInpXSecFile;           
int             gOptNev;              
int             gOptPdg = -1; //all flavors                   
double          gOptEmono = -1; //no monoenergetic
double          gOptEmin = 1e2;
double          gOptEmax = 1e10;
double          gOptCthmono = -2; //no fixed angle                   
double          gOptCthmin = -1.;
double          gOptCthmax =  1.;
double          gOptAlpha = 1; //flat spectrum in log10e
double          gOptDepth = 0.; //Detector position originally at surface
double          gOptRadius = 0.; //Detector Radius (radius of a cylindrical shape) originally a point
double          gOptHeight = 0.; //Detector height (height of the cylinder centerd at gOptDepth) originally a point
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
  IncomingFlux * flx_driver = new IncomingFlux(gOptPdg, gOptAlpha, gOptCthmin, gOptCthmax, gOptCthmono, gOptEmin, gOptEmax, gOptEmono,gOptDepth,gOptRadius,gOptHeight);

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

  LOG("ComputeAttenuation", pDEBUG) << "Initializing TAUSIC...";
  string tausicpath = gSystem->Getenv("TAUSIC");
  string cp_command = "cp "+tausicpath+"/tausic*.dat .";
  system(cp_command.c_str());

  initialize_tausic_sr_(&TAUIN);

  LOG("ComputeAttenuation", pDEBUG) << "Initializing Random...";
  RandomGen * rnd = RandomGen::Instance();

  LOG("ComputeAttenuation", pINFO) << "Creating output name";
  string OutEvFile = gOptOutDir + "/NuEarthProp_" + RunOpt::Instance()->Tune()->Name() + Form("_nu%d",gOptPdg);
  if (gOptCthmono!=-2.) OutEvFile += Form("_cth%g",gOptCthmono);
  else                  OutEvFile += Form("_cth%g-%g",gOptCthmin,gOptCthmax);
  if (gOptEmono!=-1.)   OutEvFile += Form("_e%g",gOptEmono);
  else                  OutEvFile += Form("_e%g-%g",gOptEmin,gOptEmax);
  OutEvFile += Form("_s%d",gOptRanSeed);
  OutEvFile += ".root";

  LOG("ComputeAttenuation", pNOTICE) << "@@ Output file name: " << OutEvFile;

  // create file so tree is saved inside
  TFile * outfile = new TFile(OutEvFile.c_str(),"RECREATE");
  
  int NuIn_Pdg,NuOut_Pdg;
  double NuIn_X4[4],NuOut_X4[4];
  double NuIn_P4[4],NuOut_P4[4];
  int NInt,NIntCC,NIntNC;
  int SctID[NMAXINT];
  int IntID[NMAXINT];
  int Tgt[NMAXINT];
  int Nu[NMAXINT];

  TTree *influx = new TTree("InFlux","InFlux");
  influx->Branch("Nu_In", &NuIn_Pdg, "Nu_In/I"    );       
  influx->Branch("X4_In",  NuIn_X4,  "X4_In[4]/D" );       
  influx->Branch("P4_In",  NuIn_P4,  "P4_In[4]/D" );       

  TTree *outflux = new TTree("OutFlux","OutFlux");
  outflux->Branch("Nu_In",  &NuIn_Pdg,  "Nu_In/I"       );       
  outflux->Branch("X4_In",  NuIn_X4,    "X4_In[4]/D"    );       
  outflux->Branch("P4_In",  NuIn_P4,    "P4_In[4]/D"    );       
  outflux->Branch("Nu_Out", &NuOut_Pdg, "Nu_Out/I"      );       
  outflux->Branch("X4_Out", NuOut_X4,   "X4_Out[4]/D"   );       
  outflux->Branch("P4_Out", NuOut_P4,   "P4_Out[4]/D"   );       
  outflux->Branch("NInt",   &NInt,      "NInt/I"        );       
  outflux->Branch("NIntCC", &NIntCC,    "NIntCC/I"      );       
  outflux->Branch("NIntNC", &NIntNC,    "NIntNC/I"      );       
  outflux->Branch("SctID",  SctID,      "SctID[NInt]/I" );       
  outflux->Branch("IntID",  IntID,      "IntID[NInt]/I" );       
  outflux->Branch("Tgt",    Tgt,        "Tgt[NInt]/I"   );        
  outflux->Branch("Nu",     Nu,         "Nu[NInt]/I"    );       

  LOG("ComputeAttenuation", pDEBUG) << "Event loop...";

  GHepParticle * Nu_In = new GHepParticle();
  GHepParticle * Nu_Out = new GHepParticle();
  vector<GHepParticle> SecNu;

  int NNu = 0;
  NIntCC = 0;
  NIntNC = 0;
  NInt = 0;
  while ( NNu<gOptNev ) {

    if ( SecNu.size()>0 ) {
      LOG("ComputeAttenuation", pDEBUG) << "@@SecNu      = "   << SecNu[0].Pdg() << " , E = " << SecNu[0].E();
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << SecNu[0].Vx() << " , " << SecNu[0].Vy() << " , " << SecNu[0].Vz() << " , " << SecNu[0].Vt() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << SecNu[0].Px()/SecNu[0].E() << " , " << SecNu[0].Py()/SecNu[0].E() << " , " << SecNu[0].Pz()/SecNu[0].E() << " ]";
      flx_driver->InitNeutrino(SecNu[0]);
      SecNu.erase(SecNu.begin()); //remove first entry because it is just process
    }

    EventRecord * event = trj_driver->GenerateEvent();

    if (NInt==0) {
      if ( NNu%100==0 ) LOG("ComputeAttenuation", pINFO) << "Event " << NNu << " out of " << gOptNev;
      Nu_In = flx_driver->GetNeutrino();
      FillNeutrino(Nu_In,NuIn_Pdg,NuIn_P4[0],NuIn_P4[1],NuIn_P4[2],NuIn_P4[3],NuIn_X4[0],NuIn_X4[1],NuIn_X4[2],NuIn_X4[3]);
      influx->Fill();
    }

    LOG("ComputeAttenuation", pDEBUG) << "Shooting Neutrino: " << NNu << "( " << NInt << " ) ----->" ;
    LOG("ComputeAttenuation", pDEBUG) << "@@Flux = " << flx_driver->GetNeutrino()->Pdg() << " , E =" << flx_driver->GetNeutrino()->E();
    LOG("ComputeAttenuation", pDEBUG) << "Postion   = [ " << flx_driver->GetNeutrino()->Vx() << " , " << flx_driver->GetNeutrino()->Vy() << " , " << flx_driver->GetNeutrino()->Vz() << " , " << flx_driver->GetNeutrino()->Vt() << " ] ";
    LOG("ComputeAttenuation", pDEBUG) << "Direction = [ " << flx_driver->GetNeutrino()->Px()/flx_driver->GetNeutrino()->E() << " , " << flx_driver->GetNeutrino()->Py()/flx_driver->GetNeutrino()->E() << " , " << flx_driver->GetNeutrino()->Pz()/flx_driver->GetNeutrino()->E() << " ]";

    if (event) {

      SctID[NInt] = event->Summary()->ProcInfoPtr()->ScatteringTypeId();
      IntID[NInt] = event->Summary()->ProcInfoPtr()->InteractionTypeId();
      Tgt[NInt]   = ( event->TargetNucleus() ) ? event->TargetNucleus()->Pdg() : event->HitNucleon()->Pdg();
      Nu[NInt]    = event->Probe()->Pdg();

      const TLorentzVector & X4  = *event->Vertex();
      const TLorentzVector & P4  = *event->Probe()->P4();

      // print-out
      LOG("ComputeAttenuation", pDEBUG) << "Neutrino interaction!!!";
      LOG("ComputeAttenuation", pDEBUG) << "@@Interact   = "   << SctID[NInt] << "   " << IntID[NInt];
      LOG("ComputeAttenuation", pDEBUG) << "@@Target     = "   << Tgt[NInt];
      LOG("ComputeAttenuation", pDEBUG) << "@@Probe      = "   << Nu[NInt] << " , E =" << P4.E();
      LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << X4.X() << " , " << X4.Y() << " , " << X4.Z() << " , " << X4.T() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << P4.Px()/P4.E() << " , " << P4.Py()/P4.E() << " , " << P4.Pz()/P4.E() << " ]";
      LOG("ComputeAttenuation", pDEBUG) << *event;
      
      if (IntID[NInt]==2) NIntCC++;
      else                NIntNC++;
      NInt++;

      GHepParticle * p = 0;
      TObjArrayIter piter(event);

      while ( (p=(GHepParticle *)piter.Next()) ) {
        if ( pdg::IsNeutrino(TMath::Abs(p->Pdg())) && p->Status() && p->E()>gOptEmin ) {
          p->SetEnergy( TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()) ); //not using E to avoid problems with neutrinos from pythia not on shell
          p->SetPosition( X4.X()+p->Vx()*1e-15, X4.Y()+p->Vy()*1e-15, X4.Z()+p->Vz()*1e-15, X4.T()+p->Vt() ); //position -> passing from fm [genie] to m
          SecNu.push_back(*p);
        }
        else if ( pdg::IsTau(TMath::Abs(p->Pdg())) && p->Status() ) {
          p->SetEnergy( TMath::Sqrt(p->Px()*p->Px()+p->Py()*p->Py()+p->Pz()*p->Pz()+mtau*mtau) ); //use this energy to avoid energy conservation in tauola
          p->SetPosition( X4.X()+p->Vx()*1e-15, X4.Y()+p->Vy()*1e-15, X4.Z()+p->Vz()*1e-15, X4.T()+p->Vt() ); //position -> passing from fm [genie] to m
          vector<GHepParticle> Prod = DecayTau(p);
          for (int i=0; i<Prod.size(); i++) SecNu.push_back(Prod[i]);
        }
      }

      delete event;

      if ( SecNu.size()>0 ) continue;

      LOG("ComputeAttenuation", pDEBUG) << " ----> Absorbed by the Earth";
      Nu_Out->SetPdgCode(0);
      Nu_Out->SetMomentum(0,0,0,0);
      Nu_Out->SetPosition(0,0,0,0);

    }
    else {
      LOG("ComputeAttenuation", pDEBUG) << " ----> Goodbye Earth!!!";
      Nu_Out = flx_driver->GetNeutrino();
    }

    FillNeutrino(Nu_Out,NuOut_Pdg,NuOut_P4[0],NuOut_P4[1],NuOut_P4[2],NuOut_P4[3],NuOut_X4[0],NuOut_X4[1],NuOut_X4[2],NuOut_X4[3]);
    outflux->Fill();

    if ( SecNu.size()==0 ) {
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

  system("rm geometry-Earth.root");
  system("rm tausic*.dat");

  return 0;
}

//**************************************************************************
//**************************************************************************
//THIS FUNCTION DECAY TAUS
//**************************************************************************
//**************************************************************************

vector<GHepParticle> DecayTau(GHepParticle * tau) {

  int pdgi = tau->Pdg();

  double vxi = tau->Vx()*1e-13; //from fm(genie) to cm(tausic)
  double vyi = tau->Vy()*1e-13;
  double vzi = tau->Vz()*1e-13;
  double ti  = tau->Vt()*1e9; //from s(genie) to ns(tausic)

  double momi = tau->P4()->P();
  double dxi  = tau->Px()/momi;
  double dyi  = tau->Py()/momi;
  double dzi  = tau->Pz()/momi;
  double ei   = TMath::Sqrt( momi*momi + mtau*mtau );
  
  LOG("ComputeAttenuation", pDEBUG) << "Before decay: " << pdgi << ", E = " << ei;
  LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << vxi << " , " << vyi << " , " << vzi << " , " << ti << " ]";
  LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";

  int idec; //decay flag (0=not decay // >0=decay)
  double vxf,vyf,vyf,tf; //position after propagation in meters
  double dxf,dyf,dyf,ef; //direction after propagation

  if (gOptEnableEnergyLoss) {

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

    double tauti,tautf; //dummy only used when ITFLAG!=0
    double depthf; //total distance travel after propagation

    tau_transport_sr_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&RHOR,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);

  }
  else {

    idec = 1;
    dxf = dxi;
    dyf = dyi;
    dzf = dzi;
    ef  = ef;

    // position based on lifetime
    RandomGen * rnd = RandomGen::Instance();
    double d_r = (gOptEnableDecayLength) ? -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime : 0;
    vxf = vxi + d_r*dxi*momi/mtau;
    vyf = vyi + d_r*dyi*momi/mtau;
    vzf = vzi + d_r*dzi*momi/mtau;
    tf  =  ti + d_r*ei/mtau/cSpeed;

  }

  LOG("ComputeAttenuation", pDEBUG) << "Before decay: " << pdgf << ", E = " << ef;
  LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << vxf << " , " << vyf << " , " << vzf << " , " << tf << " ]";
  LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";


  vector<GHepParticle> Prod;

  if (idec>0) {

    double momf = TMath::Sqrt( ef*ef - mtau*mtau );

    //decay tau
    TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
    TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( pdgi, 1, dxf*momf, dyf*momf, dzf*momf, ef, mtau, -1, -1, -1, -1 );
    Tauola_evt->addParticle(Tauola_tau);
    double pol = 1; //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped
    Tauola::decayOne(Tauola_tau,true,0.,0.,pol);

    for ( int sec=1; sec<Tauola_evt->getParticleCount(); sec++ ) {
      int spdg   = Tauola_evt->getParticle(sec)->getPdgID();
      double se  = Tauola_evt->getParticle(sec)->getE();
      if ( pdg::IsNeutrino(TMath::Abs(pdg)) && e>gOptEmin ) {
        double spx = Tauola_evt->getParticle(sec)->getPx();
        double spy = Tauola_evt->getParticle(sec)->getPy();
        double spz = Tauola_evt->getParticle(sec)->getPz();
        LOG("ComputeAttenuation", pDEBUG) << "Product: " << spdg << ", E = " << se;
        LOG("ComputeAttenuation", pDEBUG) << "  Position   = [ " << vxf << " , " << vyf << " , " << vzf << " , " << tf << " ]";
        LOG("ComputeAttenuation", pDEBUG) << "  Direction  = [ " << spx/se << " , " << spy/se << " , " << spz/se << " ]";
        Prod.push_back(GHepParticle(spdg,kIStUndefined,-1,-1,-1,-1,spx,spy,spz,se,vxf/100.,vyf/100.,vzf/100.,tf/1e9));
      }
    }

    delete Tauola_evt;

  }
  else {

      LOG("ComputeAttenuation", pWARN) << "Tau did not decay!!!";
      LOG("ComputeAttenuation", pWARN) << "  Energyi     = " << ei;
      LOG("ComputeAttenuation", pWARN) << "  Positioni   = [ " << vxi << " , " << vyi << " , " << vzi << " , " << ti << " ]";
      LOG("ComputeAttenuation", pWARN) << "  Directioni  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";
      LOG("ComputeAttenuation", pWARN) << "  Energyf     = " << ef;
      LOG("ComputeAttenuation", pWARN) << "  Positionf   = [ " << vxf << " , " << vyf << " , " << vzf << " , " << tf << " ]";
      LOG("ComputeAttenuation", pWARN) << "  Directionf  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";    

  }

  return Prod;

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

  const int NLAYERS = 10;
  TString fComp[NLAYERS] = { "Core", "Core", "Mantle", "Mantle", "Mantle", "Mantle", "Mantle", "Mantle", "Rock", "Ice"       };
  double fR[NLAYERS]     = { 1221.5,  3480.,    5701.,    5771.,    5971.,    6151.,   6346.6,    6356.,  6368.,  fREarth_km }; //km

  double fEarthCoeff[NLAYERS][4] = {
    {13.0885,0.,-8.8381,0.},
    {12.5815,-1.2638,-3.6426,-5.5281} ,
    {7.9565,-6.4761,5.5283,-3.0807},
    {5.3197,-1.4836,0.,0.},
    {11.2494,-8.0298,0.,0.},
    {7.1089,-3.8045,0.,0.},
    {2.691,0.6924,0.,0.},
    {2.9,0.,0.,0.},
    {2.6,0.,0.,0.},
    {0.9168,0.,0.,0.}
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
    
    if( r1>=fR[NLAYERS-4] ) break;
    
  }

  for(int iiLayer=NLAYERS-3; iiLayer<NLAYERS; iiLayer++){ // constant density
    Layer layer;
    layer.r1          = fR[iiLayer-1];
    layer.r2          = fR[iiLayer];
    layer.rho         = fEarthCoeff[iiLayer][0] + fEarthCoeff[iiLayer][1]*layer.r2/fREarth_km + fEarthCoeff[iiLayer][2]*pow(layer.r2/fREarth_km,2) + fEarthCoeff[iiLayer][3]*pow(layer.r2/fREarth_km,3);
    layer.Composition = fComp[iiLayer];    
    VecLayers.push_back(layer);
  }

  map<int,double> RockComp, MantleComp, CoreComp, IceComp;
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
  IceComp     .insert(map<int, double>::value_type(1000080160,0.8881)    );
  IceComp     .insert(map<int, double>::value_type(1000010010,0.1119)    );

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
    if ( VecLayers[iiLayer].Composition=="Ice" )     {
      LayerMix[iiLayer] = new TGeoMixture( name, IceComp.size(), VecLayers[iiLayer].rho ); 
      iter = IceComp.begin();
      for( ; iter != IceComp.end(); ++iter )     LayerMix[iiLayer]->AddElement( pdg::IonPdgCodeToA(iter->first), pdg::IonPdgCodeToZ(iter->first), iter->second );                  
    }    
    else if ( VecLayers[iiLayer].Composition=="Rock" )     {
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
      if(opt.compare("--depth")==0){
        i++; 
        LOG("ComputeAttenuation", pINFO) << "Reading detector depth";
        gOptDepth = atof(argv[i]);  
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
     << "\n\t Depth        = " << gOptDepth << " m"
     << "\n\t Radius       = " << gOptRadius << " m"
     << "\n\t Height       = " << gOptHeight << " m"
     << "\n\t EnergyLoss   = " << gOptEnableEnergyLoss
     << "\n\t DecayLength  = " << gOptEnableDecayLength
     << "\n\n";


}



