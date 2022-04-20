#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include <TString.h> 
#include <TSystem.h> 

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

#include <map>

using namespace std;
using namespace genie;


const double fPREM_ref = 6371.0;

struct XmlLayer{ 
  TString name; 
  double radius; 
  double density[4]; 
  map<int,double> elements;
};

struct Layer{ 
  TString name; 
  double r1,r2; 
  double rho; 
  map<int,double> elements;
};


void GetCommandLineArgs (int argc, char ** argv);

string pathxml;
string pathroot;

int main (int argc, char** argv) 
{

  // Parse command line arguments
  LOG("BuildEarth", pDEBUG) << "Reading options...";
  GetCommandLineArgs(argc,argv);

  LOG("BuildEarth", pDEBUG) << pathxml << "  " << pathroot;

  std::vector<XmlLayer> xmllayers;

  int auxpdg = 0;

  xmlTextReaderPtr reader = xmlNewTextReaderFilename(pathxml.c_str());
  while(xmlTextReaderRead(reader)==1){

    xmlChar * name  = xmlTextReaderName     (reader);
    xmlChar * value = xmlTextReaderValue    (reader);
    int       type  = xmlTextReaderNodeType (reader);
    int       depth = xmlTextReaderDepth    (reader);
    
    if(depth==0 && type==1){
      if(xmlStrcmp(name, (const xmlChar *) "earth_prem")){
      LOG("BuildEarth", pFATAL) << " invalid root element in filename " << pathxml <<"; exiting... ";
      exit(1);
     }
    }

    if(!xmlStrcmp(name, (const xmlChar*) "param_set" ) && type==1){

      XmlLayer auxlayer;

      xmlChar * xname   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
      auxlayer.name = utils::str::TrimSpaces((const char *)xname);

      xmlChar * xradius = xmlTextReaderGetAttribute(reader,(const xmlChar*)"radius");
      auxlayer.radius = stof(utils::str::TrimSpaces((const char *)xradius));

      xmlChar * xdensity = xmlTextReaderGetAttribute(reader,(const xmlChar*)"density");
      std::vector<string> sdensity = utils::str::Split(utils::str::TrimSpaces((const char *)xdensity), ",");
      assert(sdensity.size()==4);
      for ( int i=0 ; i<4; i++ ) auxlayer.density[i] = stof(sdensity[i]);

      xmllayers.push_back(auxlayer);
    
    }

    if(!xmlStrcmp(name, (const xmlChar *) "param" ) && type==1) {
      xmlChar * xelem   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"element");
      auxpdg = stoi(utils::str::TrimSpaces((const char *)xelem));
    }

    if(depth==3) xmllayers.back().elements.insert({auxpdg,atof((const char *)value)});

  }


  double fREarth_km = xmllayers.back().radius;

  LOG("BuildEarth", pDEBUG) << fREarth_km;

  std::vector<Layer> layers;
  
  double rmean = 0.;
  double step  = 100.;
  double r1    = 0;
  double r2    = step;  
  bool border  = false;  
  int iLayer   = 0;
  int iStep    = 0;
    
  while(1){       

    if( r2>xmllayers[iLayer].radius || (xmllayers[iLayer].density[1]==0 && xmllayers[iLayer].density[2]==0 && xmllayers[iLayer].density[3]==0) ) {
      r2     = xmllayers[iLayer].radius; 
      border = true;
    }

    rmean             = ( r2 + r1 ) / 2.;

    Layer auxlayer;
    auxlayer.name        = xmllayers[iLayer].name;
    auxlayer.r1          = r1;
    auxlayer.r2          = r2;
    auxlayer.rho         = xmllayers[iLayer].density[0] + xmllayers[iLayer].density[1]*rmean/fPREM_ref + xmllayers[iLayer].density[2]*pow(rmean/fPREM_ref,2) + xmllayers[iLayer].density[3]*pow(rmean/fPREM_ref,3);
    auxlayer.elements    = xmllayers[iLayer].elements;
    
    layers.push_back(auxlayer);
    
    if( border ) {
      r1     = r2;
      border = false;
      iLayer++;
    }
    else r1 = (iStep+1)*step;

    r2 = (iStep+2)*step;
    iStep++;
 
    if( r1>=fREarth_km ) break;

  }

  // Define media and volumes
  TGeoManager * GeoManager =  new TGeoManager("VolGenGeo", "generation volume geometry");
  TGeoTranslation * Trans = new TGeoTranslation(0.,0.,0.);

  // vacuum
  TGeoMaterial * MatVacuum = new TGeoMaterial("MatVacuum");
  MatVacuum->SetA(0.);
  MatVacuum->SetZ(0.);
  MatVacuum->SetDensity(0.);
  TGeoMedium * Vacuum = new TGeoMedium("Vacuum", 0, MatVacuum);
  TGeoVolume *TopVolume = GeoManager->MakeSphere("TopVolume",Vacuum,0.,10.*layers.back().r2*1.E3);
  GeoManager->SetTopVolume(TopVolume);

  int counter = 0;  
  for ( auto layer : layers ) {
    LOG("BuildEarth", pDEBUG) << layer.name << ", r1=" << layer.r1 << ", r2=" << layer.r2 << ", rho=" << layer.rho;
    TGeoMixture * LayerMix = new TGeoMixture( layer.name, layer.elements.size(), layer.rho ); 
    for( auto elem : layer.elements ) {    
      LOG("BuildEarth", pDEBUG) << elem.first << "  " << elem.second;
      LayerMix->AddElement( pdg::IonPdgCodeToA(elem.first), pdg::IonPdgCodeToZ(elem.first), elem.second );                  
    }
    TGeoMedium * LayerMedium = new TGeoMedium( layer.name, counter+1, LayerMix );   
    TGeoVolume * Sphere = GeoManager->MakeSphere( layer.name, LayerMedium, layer.r1*1.E3, layer.r2*1.E3 ); 
    TopVolume->AddNode( Sphere, counter+1, Trans );
    counter++;
  }

  GeoManager->CloseGeometry();
  GeoManager->Export(pathroot.c_str());

  return 0;

}

void GetCommandLineArgs(int argc, char ** argv)
{

  size_t apos;
  string opt;

  for (int i = 1; i < argc; i++) {
    string arg(argv[i]);
    if (arg.compare(0,1,"-")==0){
      apos=arg.find(" ");
      opt=arg.substr(0,apos);

      if(opt.compare("-xml")==0){   
        i++;
        LOG("BuildEarth", pDEBUG) << "Reading xml file";
        pathxml = argv[i];
        if ( gSystem->AccessPathName(pathxml.c_str()) ) {
          LOG("BuildEarth", pFATAL) << "Xml file not found: " << pathxml << "; exit";
          exit(1);
        }
      }          
      if(opt.compare("-root")==0){   
        i++;
        LOG("BuildEarth", pDEBUG) << "Reading output root file";
        pathroot = argv[i];
      }          
    }
  }

}



