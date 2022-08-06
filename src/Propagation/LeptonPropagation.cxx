#include "LeptonPropagation.h"


std::map<int, int> proposal_translation = {
  {1000000002,11},
  {1000000003,11},
  {1000000004,11},
  {1000000005,211},
  {1000000008,-1},
};

LeptonPropagation::LeptonPropagation(int pdg, string proposaltable, double ecut, double vcut, int seed, ROOTGeomAnalyzer * gd, vector<string> skiplist) {

  geom_driver = gd;

  PROPOSAL::RandomGenerator::Get().SetSeed(seed);

  ConfigProposal(pdg,proposaltable,ecut,vcut,skiplist);

}

std::vector<PROPOSAL::Components::Component> LeptonPropagation::GetComponent(map<int,double> composition) {

  NaturalIsotopes * iso = NaturalIsotopes::Instance();

  double MassFracMax=0.;
  double MoliFracMass=1.;
  for ( auto c : composition ) {
    int    PdgCode  = c.first;
    int    Z        = pdg::IonPdgCodeToZ(PdgCode);
    double MassFrac = c.second;
    double AtomMass = iso->ElementDataPdg(Z,PdgCode)->AtomicMass();
    if ( MassFrac>MassFracMax ) {
      MassFracMax  = MassFrac;
      MoliFracMass = MassFrac/AtomMass;
    }
  }

  std::vector<PROPOSAL::Components::Component> component;
  for ( auto c : composition ) {
    int    PdgCode  = c.first;
    int    Z        = pdg::IonPdgCodeToZ(PdgCode);
    double MassFrac = c.second;
    double AtomMass = iso->ElementDataPdg(Z,PdgCode)->AtomicMass();
    double AtomInMolecule = MassFrac/AtomMass/MoliFracMass;     

    component.push_back(PROPOSAL::Components::Component(to_string(PdgCode).c_str(),Z,AtomMass,AtomInMolecule));
  
  }

  return component;

}


void LeptonPropagation::ConfigProposal(int pdg, string proposaltable, double ecut, double vcut, vector<string> skiplist ) {

  std::vector<shared_ptr<const PROPOSAL::Geometry>> vgeolayers;
  std::vector<shared_ptr<const PROPOSAL::Medium>> vmedlayers;

  //numbers extracted from proposal medium definitions
  map<string,Ionisation_Constants> ion_const;
  ion_const["Ice"]    = {  75.0, -3.50170, 0.09116, 3.4773,  0.2400, 2.8004, 0.00 };
  ion_const["Mantle"] = { 136.4, -3.77380, 0.08301, 3.4120,  0.0492, 3.0549, 0.00 };
  ion_const["Core"]   = { 286.0, -4.29110, 0.14680, 2.9632, -0.0012, 3.1531, 0.12 };
  ion_const["Atmos"]  = {  85.7, -10.5961, 0.10914, 3.3994,  1.7418, 4.2759, 0.00 };

  TLorentzVector p4(0,0,1,1);
  TLorentzVector x4(0,0,0,0);
  std::vector< pair<double, const TGeoMaterial*> > MatLengths = geom_driver->ComputeMatLengths(x4,p4);

  double lin  = 0;
  double lout = 0;
  for ( auto sitr : MatLengths ) {

    const  TGeoMaterial* mat = sitr.second;
    double mat_rho = mat->GetDensity(); //gr/cm3

    lin = lout;
    lout += sitr.first; //m
    string mat_name = mat->GetName();

    auto glayer = make_shared<PROPOSAL::Sphere>(PROPOSAL::Vector3D(0,0,0), lout, lin);
    glayer->SetHierarchy(1);    

    bool skip = false;
    for ( auto s : skiplist ) {
      if ( s==mat->GetName() ) { 
        LOG("LeptonPropagation", pWARN) << "Skipping material for PROPOSAL: " << mat->GetName() << " [ " << lin << " , " << lout << " ]";
        skip = true; break; 
      }
    }

    if ( skip ) continue;

    map<int,double> composition;
    if (mat->IsMixture()) {
      const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
      for (int i = 0; i < mixt->GetNelements(); i++) composition[geom_driver->GetTargetPdgCode(mixt, i)] = mixt->GetWmixt()[i];
    }
    else composition[geom_driver->GetTargetPdgCode(mat)] = 1.;

    string lname = "";
    if      (mat_name.find("Ice")    != std::string::npos) lname="Ice";
    else if (mat_name.find("Mantle") != std::string::npos) lname="Mantle";
    else if (mat_name.find("Core")   != std::string::npos) lname="Core";
    else if (mat_name.find("Atmos")  != std::string::npos) lname="Atmos";
    
    auto mlayer = make_shared<const PROPOSAL::Medium>(mat_name, 1, ion_const[lname].I, ion_const[lname].C, ion_const[lname].a, ion_const[lname].m, ion_const[lname].X0, ion_const[lname].X1, ion_const[lname].d0, mat_rho, GetComponent(composition));
    vmedlayers.push_back(mlayer);     
    vgeolayers.push_back(glayer);


  }

  auto geoWorld = make_shared<PROPOSAL::Sphere>(PROPOSAL::Vector3D(0,0,0), lout, 0 );

  unique_ptr<PROPOSAL::Sector::Definition> def_global(new PROPOSAL::Sector::Definition());

  PROPOSAL::Sector::Definition sec_def=*def_global;
  
  LOG("LeptonPropagation", pWARN) << "Ecut: " << ecut << "  Vcut: " << vcut;
  if (ecut>0) sec_def.cut_settings.SetEcut(ecut*1e3);
  else        sec_def.cut_settings.SetEcut(-1);
  sec_def.cut_settings.SetVcut(vcut);
  // sec_def.do_continuous_randomization = false;    
  sec_def.do_continuous_energy_loss_output = false;    

  //matching danton  
  // sec_def.scattering_model = PROPOSAL::ScatteringFactory::Moliere;
  // sec_def.utility_def.photo_def.shadow = PROPOSAL::PhotonuclearFactory::ShadowDuttaRenoSarcevicSeckel;
  // sec_def.utility_def.photo_def.parametrization = PROPOSAL::PhotonuclearFactory::AbramowiczLevinLevyMaor97;
  // sec_def.utility_def.brems_def.parametrization = PROPOSAL::BremsstrahlungFactory::SandrockSoedingreksoRhode;
  // sec_def.utility_def.brems_def.lpm_effect = false;
  // sec_def.utility_def.epair_def.parametrization = PROPOSAL::EpairProductionFactory::SandrockSoedingreksoRhode;
  // sec_def.utility_def.epair_def.lpm_effect = false;

  std::vector<PROPOSAL::Sector::Definition> sectors;
  for (unsigned int i=0; i<vgeolayers.size(); i++ ) {
    sec_def.SetMedium(vmedlayers[i]);
    sec_def.SetGeometry(vgeolayers[i]);
    // cout << sec_def << endl;
    sectors.push_back(sec_def);
  }

  PROPOSAL::InterpolationDef interpolation_def;
  interpolation_def.path_to_tables          = proposaltable;
  interpolation_def.path_to_tables_readonly = proposaltable;
  interpolation_def.nodes_cross_section = 200;

  
  LOG("LeptonPropagation", pNOTICE) << "Loading Proposal tables... " << pdg << " -> " << proposaltable;
  if      (pdg==13) Proposal = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(),  sectors, geoWorld, interpolation_def );
  else if (pdg==15) Proposal = new PROPOSAL::Propagator(PROPOSAL::TauMinusDef::Get(), sectors, geoWorld, interpolation_def );

  return;

}

double LeptonPropagation::ComputeDepth(GHepParticle * p) { //length(m)

  std::vector< std::pair<double, const TGeoMaterial*> > MatLengths = geom_driver->ComputeMatLengths(*p->X4(),*p->P4());

  if (MatLengths[0].second->GetDensity()==0) {
    LOG("LeptonPropagation", pFATAL) << "Lepton can not start in vacum (neutrinos do not interact!)";
    exit(1);
  }

  double totlength = 0.;
  for ( auto sitr : MatLengths ) {
    if (sitr.second->GetDensity()>0) totlength += sitr.first; //only propagate until no vacum region
  }
  LOG("LeptonPropagation", pDEBUG) << "Sum ---> Length = " << totlength << " m";

  return totlength;

}


void LeptonPropagation::Step(GHepParticle * lepton, double length, double minenergy) //m and GeV
{ 

  int pdgi = lepton->Pdg();

  double vxi = lepton->Vx()*1e2; //from m(geom) to (proposal)
  double vyi = lepton->Vy()*1e2;
  double vzi = lepton->Vz()*1e2;
  double ti  = lepton->Vt();

  double momi = lepton->P4()->P();
  double dxi  = lepton->Px()/momi;
  double dyi  = lepton->Py()/momi;
  double dzi  = lepton->Pz()/momi;
  double ei   = TMath::Sqrt( momi*momi + mass*mass );
  
  LOG("LeptonPropagation", pDEBUG) << "Before step: " << pdgi << ", E = " << ei << " GeV";
  LOG("LeptonPropagation", pDEBUG) << "  Position   = [ " << vxi << " cm, " << vyi << " cm, " << vzi << " cm, " << ti << " ns ]";
  LOG("LeptonPropagation", pDEBUG) << "  Direction  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";

  PROPOSAL::Vector3D direction(dxi,dyi,dzi);
  direction.CalculateSphericalCoordinates();

  PROPOSAL::DynamicData particle(TMath::Abs(pdgi));

  particle.SetEnergy(ei*1E3); // [MeV]
  particle.SetPropagatedDistance(0);
  particle.SetPosition(PROPOSAL::Vector3D(vxi,vyi,vzi));
  particle.SetDirection(direction);
  particle.SetTime(ti);

  PROPOSAL::Secondaries sec = Proposal->Propagate(particle,length*1e2,minenergy*1e3); //cm and MeV

  auto particles = sec.GetModifyableSecondaries();

  int nparticles = particles.size();

  GHepStatus_t idec = kIStUndefined; //decay flag
  double e_dec = 0.;
  for (int i=nparticles-1; i>=0; i--) {
    if ( particles[i].GetType()>1000000000 ) break;
    else {
      idec = kIStDecayedState;
      e_dec += particles[i].GetEnergy();
    }
  }

  double ef;
  if ( idec==kIStDecayedState ) ef = e_dec/1e3; //from MeV to GeV
  else                          ef = sec.GetEnergy().back()/1e3; //from MeV to GeV

  double vxf = particles.back().GetPosition().GetX();
  double vyf = particles.back().GetPosition().GetY();
  double vzf = particles.back().GetPosition().GetZ();
  double dxf = particles.back().GetDirection().GetX();
  double dyf = particles.back().GetDirection().GetY();
  double dzf = particles.back().GetDirection().GetZ();
  double tf  = particles.back().GetTime();

  // ef  = ei*0.99;
  // idec = (ef<1e2) ? kIStDecayedState : kIStUndefined;
  // vxf = vxi+length*1e2*dxi;
  // vyf = vyi+length*1e2*dyi;
  // vzf = vzi+length*1e2*dzi;
  // dxf = dxi;
  // dyf = dyi;
  // dzf = dzi;
  // tf  = ti;

  LOG("LeptonPropagation", pDEBUG) << "After step: " << pdgi << ", E = " << ef << " GeV" << ", idec = " << idec;
  LOG("LeptonPropagation", pDEBUG) << "  length     = " << sqrt(pow(vxi-vxf,2)+pow(vyi-vyf,2)+pow(vzi-vzf,2)) << " cm";
  LOG("LeptonPropagation", pDEBUG) << "  Position   = [ " << vxf << " cm, " << vyf << " cm, " << vzf << " cm, " << tf << " ns ]";
  LOG("LeptonPropagation", pDEBUG) << "  Direction  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";

  vxf *= 1e-2; //from cm(proposal) to m(geom)
  vyf *= 1e-2;
  vzf *= 1e-2;

  double momf = TMath::Sqrt( ef*ef - mass*mass );

  lepton->SetStatus(idec);
  lepton->SetMomentum(dxf*momf,dyf*momf,dzf*momf,ef);
  lepton->SetPosition(vxf,vyf,vzf,tf);

}


std::vector<GHepParticle> LeptonPropagation::StepShowers(GHepParticle * lepton, double length, double minenergy) //m and GeV
{ 

  int pdgi = lepton->Pdg();

  double vxi = lepton->Vx()*1e2; //from m(geom) to (proposal)
  double vyi = lepton->Vy()*1e2;
  double vzi = lepton->Vz()*1e2;
  double ti  = lepton->Vt();

  double momi = lepton->P4()->P();
  double dxi  = lepton->Px()/momi;
  double dyi  = lepton->Py()/momi;
  double dzi  = lepton->Pz()/momi;
  double ei   = TMath::Sqrt( momi*momi + mass*mass );
  
  LOG("LeptonPropagation", pDEBUG) << "Before step: " << pdgi << ", E = " << ei << " GeV";
  LOG("LeptonPropagation", pDEBUG) << "  Position   = [ " << vxi << " cm, " << vyi << " cm, " << vzi << " cm, " << ti << " s ]";
  LOG("LeptonPropagation", pDEBUG) << "  Direction  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";

  PROPOSAL::Vector3D direction(dxi,dyi,dzi);
  direction.CalculateSphericalCoordinates();

  PROPOSAL::DynamicData particle(TMath::Abs(pdgi));

  particle.SetEnergy(ei*1E3); // [MeV]
  particle.SetPropagatedDistance(0);
  particle.SetPosition(PROPOSAL::Vector3D(vxi,vyi,vzi));
  particle.SetDirection(direction);
  particle.SetTime(ti);

  PROPOSAL::Secondaries sec = Proposal->Propagate(particle,length*1e2,minenergy*1e3); //cm and MeV

  auto particles = sec.GetModifyableSecondaries();

  int nparticles = particles.size();

  std::vector<GHepParticle> showers;
  for (int i=0; i<nparticles; i++) {
    double se = (particles[i].GetParentParticleEnergy()-particles[i].GetEnergy())*1e-3;

    if ( se>minenergy ) {

      int spdg = (particles[i].GetType()>1000000000) ? proposal_translation[particles[i].GetType()] : particles[i].GetType();

      if      ( pdg::IsNeutrino(TMath::Abs(spdg)) ) continue;
      else if ( spdg==-1 ) continue;
      else if ( spdg==0 ) {
        LOG("LeptonPropagation", pFATAL) << "Wrong PDG!!!";
        LOG("LeptonPropagation", pFATAL) << particles[i].GetType() << ", E = " << se << " GeV";
        exit(1);
      }
      
      double spx = particles[i].GetDirection().GetX()*se;
      double spy = particles[i].GetDirection().GetY()*se;
      double spz = particles[i].GetDirection().GetZ()*se;
      double svx = particles[i].GetPosition().GetX()*1e-2;
      double svy = particles[i].GetPosition().GetY()*1e-2;
      double svz = particles[i].GetPosition().GetZ()*1e-2;
      double svt = particles[i].GetTime();
      showers.push_back(GHepParticle(spdg,kIStUndefined,-1*TMath::Abs(particles[i].GetType()),-1,-1,-1,spx,spy,spz,se,svx,svy,svz,svt));
    }
  }

  return showers;

}


