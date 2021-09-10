
#include "TauPropagation.h"

extern "C" { 
  void initialize_tausic_sr_(int *muin);
  void initialize_tausic_sw_(int *muin);
  void tau_transport_sr_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
  void tau_transport_sw_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
}

TauPropagation::TauPropagation(string ptype, int seed, GeomAnalyzerI * gd) {

  tauproptype = ptype;

  geom_driver = gd;

  LOG("TauPropagation", pDEBUG) << "Initializing TAUOLA...";
  Tauola::initialize();
  Tauola::setSeed(seed,0,0);
  mtau         = Tauola::getTauMass();        //use this mass to avoid energy conservation warning in tauola
  d_lifetime   = Tauola::tau_lifetime/(constants::kLightSpeed/(units::millimeter/units::nanosecond)); //lifetime in ns
  polarization = 1;                           //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped
  LOG("TauPropagation", pDEBUG) << "d_lifetime (TAUOLA): " << d_lifetime;


  /*TAUSIC:
  REFERENCES:
  https://arxiv.org/pdf/hep-ph/9705408.pdf
  https://arxiv.org/pdf/hep-ph/9911493.pdf
  https://arxiv.org/pdf/0810.4635.pdf

  MODELS:
  -bremsstrahlung  ->  S.R. Kelner, R.P. Kokoulin, A.A. Petrukhin. Phys. At. Nucl. 60 (1997) 576.
  -pair production ->  R.P. Kokoulin, A.A. Petrukhin. In: Proc. 12th Intern. Cosmic Ray Conf. (Hobart, 1971), Vol. 6, p. 2436. S.R. Kelner. Phys. At. Nucl. 61 (1998) 448.
  -Inelastic       -> [BS] E.V. Bugaev, Yu.V. Shlepin. Phys. Rev. D 67 (2003) 034027.
                      [ALLM] Abramowicz, A. Levy. Preprint DESY 97-251 (1997).
  -scat. angle in bremss. and pp -> A. van Ginneken. Nucl. Istrum. and Meth. in Phys. Res. A 251 (1986) 21.
  -scat. angle in inelastic -> using diff xsec from: [BS] L.B. Bezrukov, E.V. Bugaev. Sov. J. Nucl. Phys. 32 (1980) 635; 33 (1981) 847.
                                                     [ALLM] Abramowicz, A. Levy. Preprint DESY 97-251 (1997).
  */
  if ( tauproptype=="TAUSIC-ALLM" || tauproptype=="TAUSIC-BS" ) {
    if      ( tauproptype=="TAUSIC-BS"   ) TAUIN = 1;
    else if ( tauproptype=="TAUSIC-ALLM" ) TAUIN = 2;
    TAUMODEL = 1; //mult.scatt. and scatt. due to other processes switch ON
    ITFLAG = 1; //set lifetime outside tausic
    LOG("TauPropagation", pDEBUG) << "Initializing TAUSIC...";
    string tausicpath = gSystem->Getenv("TAUSIC");
    string cp_command = "cp "+tausicpath+"/tausic*.dat .";
    system(cp_command.c_str());
    initialize_tausic_sr_(&TAUIN);      
    initialize_tausic_sw_(&TAUIN);      
    system("rm ./tausic*.dat");
  }

  else if ( tauproptype=="PROPOSAL" ) {
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
    ConfigProposal(geom_driver);
  }

  LOG("TauPropagation", pDEBUG) << "Initializing Random...";
  rnd = RandomGen::Instance();

}

std::vector<PROPOSAL::Components::Component> TauPropagation::GetComponent(map<int,double> composition) {

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


void TauPropagation::ConfigProposal(GeomAnalyzerI * gd) {

  std::vector<shared_ptr<const PROPOSAL::Geometry>> vgeolayers;
  std::vector<shared_ptr<const PROPOSAL::Medium>> vmedlayers;

  //numbers extracted from proposal medium definitions
  map<string,Ionisation_Constants> ion_const;
  ion_const["Ice"]    = {  75.0, -3.5017, 0.09116, 3.4773,  0.2400, 2.8004, 0.00 };
  ion_const["Mantle"] = { 136.4, -3.7738, 0.08301, 3.4120,  0.0492, 3.0549, 0.00 };
  ion_const["Core"]   = { 286.0, -4.2911, 0.14680, 2.9632, -0.0012, 3.1531, 0.12 };

  TLorentzVector p4(0,0,1,1);
  TLorentzVector x4(0,0,0,0);
  std::vector< pair<double, const TGeoMaterial*> > MatLengths = gd->ComputeMatLengths(x4,p4);

  double lin  = 0;
  double lout = 0;
  for ( auto sitr : MatLengths ) {

    lout += sitr.first; //m
    const  TGeoMaterial* mat = sitr.second;
    string mat_name = mat->GetName();
    double mat_rho = mat->GetDensity(); //gr/cm3

    auto glayer = make_shared<PROPOSAL::Sphere>(PROPOSAL::Vector3D(0,0,0), lout, lin);
    glayer->SetHierarchy(1);    
    vgeolayers.push_back(glayer);

    lin = lout;

    map<int,double> composition;
    if (mat->IsMixture()) {
      const TGeoMixture * mixt = dynamic_cast <const TGeoMixture*> (mat);
      for (int i = 0; i < mixt->GetNelements(); i++) composition[gd->GetTargetPdgCode(mixt, i)] = mixt->GetWmixt()[i];
    }
    else composition[gd->GetTargetPdgCode(mat)] = 1.;

    string lname = "";
    if      (mat_name.find("Ice")    != std::string::npos) lname="Ice";
    else if (mat_name.find("Mantle") != std::string::npos) lname="Mantle";
    else if (mat_name.find("Core")   != std::string::npos) lname="Core";
    
    auto mlayer = make_shared<const PROPOSAL::Medium>(mat_name, 1, ion_const[lname].I, ion_const[lname].C, ion_const[lname].a, ion_const[lname].m, ion_const[lname].X0, ion_const[lname].X1, ion_const[lname].d0, mat_rho, GetComponent(composition));
    vmedlayers.push_back(mlayer);     

  }

  auto geoWorld = make_shared<PROPOSAL::Sphere>(PROPOSAL::Vector3D(0,0,0), lout, 0 );

  unique_ptr<PROPOSAL::Sector::Definition> def_global(new PROPOSAL::Sector::Definition());

  PROPOSAL::Sector::Definition sec_def=*def_global;
  sec_def.cut_settings.SetEcut(-1);
  sec_def.cut_settings.SetVcut(0.001);

  std::vector<PROPOSAL::Sector::Definition> sectors;
  for (unsigned int i=0; i<vgeolayers.size(); i++ ) {
    sec_def.SetMedium(vmedlayers[i]);
    sec_def.SetGeometry(vgeolayers[i]);
    sectors.push_back(sec_def);
  }

  PROPOSAL::InterpolationDef interpolation_def;
  interpolation_def.path_to_tables          = string(gSystem->Getenv("NUPROPEARTH"))+"/proposal_tables";
  interpolation_def.path_to_tables_readonly = string(gSystem->Getenv("NUPROPEARTH"))+"/proposal_tables";
  interpolation_def.nodes_cross_section = 200;

  ProposalTau = new PROPOSAL::Propagator(PROPOSAL::TauMinusDef::Get(), sectors, geoWorld, interpolation_def );

  return;

}

void TauPropagation::ComputeDepth(GHepParticle * p, double &avgrho, double &totlength) { //length(cm) depth(g/cm2)

  totlength = 0.;
  avgrho    = 0.;
  std::vector< std::pair<double, const TGeoMaterial*> > MatLengths = geom_driver->ComputeMatLengths(*p->X4(),*p->P4());
  for ( auto sitr : MatLengths ) {
    double length  = sitr.first * 1e2;  //from m(geom) to cm(tau)
    double rho     = sitr.second->GetDensity();
    LOG("TauPropagation", pDEBUG) << "  Rho = " << rho << " g/cm^3 ; Length = " << length << " cm";
    totlength += length;
    avgrho    += length*rho;
  }
  avgrho = avgrho/totlength;
  LOG("TauPropagation", pDEBUG) << "Sum ---> Avg. Rho = " << avgrho << " g/cm^3 ; Length = " << totlength << " cm";

}


std::vector<GHepParticle> TauPropagation::Propagate(GHepParticle * tau) { 

  int pdgi = tau->Pdg();

  double vxi = tau->Vx()*1e2; //from m(geom) to cm(tausic)
  double vyi = tau->Vy()*1e2;
  double vzi = tau->Vz()*1e2;
  double ti  = tau->Vt()*1e9; //from s(geom) to ns(tausic)

  double momi = tau->P4()->P();
  double dxi  = tau->Px()/momi;
  double dyi  = tau->Py()/momi;
  double dzi  = tau->Pz()/momi;
  double ei   = TMath::Sqrt( momi*momi + mtau*mtau );
  
  LOG("TauPropagation", pDEBUG) << "Before decay: " << pdgi << ", E = " << ei << " GeV";
  LOG("TauPropagation", pDEBUG) << "  Position   = [ " << vxi << " cm, " << vyi << " cm, " << vzi << " cm, " << ti << " ns ]";
  LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << dxi << " , " << dyi << " , " << dzi << " ]";

  int idec = 0; //decay flag (0=not decay // 1=decay)
  double vxf,vyf,vzf,tf; //position after propagation in cm/ns
  double dxf,dyf,dzf,ef; //direction after propagation

  double avgrho,depthi; //g/cm3, cm
  ComputeDepth(tau,avgrho,depthi);

  if      ( tauproptype=="TAUSIC-ALLM" or tauproptype=="TAUSIC-BS" ) {

    double tauti = -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime; //ns
    
    while ( depthi>1 ) { //threshold in 1cm
      double tautf; //remaining lifetime [ns]
      double depthf; //total distance travel after propagation [cm]
      if   ( avgrho>2. ) tau_transport_sr_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&avgrho,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);  
      else               tau_transport_sw_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&avgrho,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);

      if ( idec==1 ) break; //tau has decayed

      double momf = TMath::Sqrt( ef*ef - mtau*mtau );
      tau->SetMomentum(dxf*momf, dyf*momf, dzf*momf, ef);
      tau->SetPosition(vxf*1e-2,vyf*1e-2,vzf*1e-2,tf*1e-9); //from cm [tausic] to m [geom]
      ComputeDepth(tau,avgrho,depthi);
      
      tauti = tautf;
      vxi = vxf; vyi = vyf; vzi = vzf; ti = tf;
      dxi = dxf; dyi = dyf; dzi = dzf; ei = ef; 
      
    }

  }
  else if ( tauproptype=="PROPOSAL" ) {

    PROPOSAL::Vector3D direction(dxi,dyi,dzi);
    direction.CalculateSphericalCoordinates();

    PROPOSAL::DynamicData particle_tau(PROPOSAL::TauMinusDef::Get().particle_type);

    particle_tau.SetEnergy(ei*1E3); // [MeV]
    particle_tau.SetPropagatedDistance(0);
    particle_tau.SetPosition(PROPOSAL::Vector3D(vxi,vyi,vzi));
    particle_tau.SetDirection(direction);
    particle_tau.SetTime(ti);

    PROPOSAL::Secondaries sec = ProposalTau->Propagate(particle_tau,depthi);

    auto particles = sec.GetSecondaries();

    int nparticles = particles.size();
    LOG("TauPropagation", pDEBUG) << "nparticles: " << nparticles;

    double e_dec = 0.;
    for (int i=nparticles-1; i>=0; i--) {
      if ( particles[i].GetType()>1000000000 ) break;
      else {
        idec = 1;
        e_dec += particles[i].GetEnergy();
      }
    }

    if ( idec==1 ) ef = e_dec/1e3; //from MeV to GeV
    else           ef = sec.GetEnergy().back()/1e3; //from MeV to GeV

    vxf = particles.back().GetPosition().GetX();
    vyf = particles.back().GetPosition().GetY();
    vzf = particles.back().GetPosition().GetZ();
    dxf = particles.back().GetDirection().GetX();
    dyf = particles.back().GetDirection().GetY();
    dzf = particles.back().GetDirection().GetZ();
    tf  = particles.back().GetTime();


  }
  else if (tauproptype=="NOELOSS") {

    double tauti = -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime; //ns

    dxf = dxi; dyf = dyi; dzf = dzi; ef = ei;
    double d_r  = tauti*(constants::kLightSpeed/(units::centimeter/units::nanosecond)); //cm
    if ( d_r*momi/mtau>depthi ) {
      idec = 0;
      vxf = vxi + depthi*dxi;
      vyf = vyi + depthi*dyi;
      vzf = vzi + depthi*dzi;
      tf  =  ti + depthi/(constants::kLightSpeed/(units::centimeter/units::nanosecond));      
    }
    else {
      idec = 1;
      vxf = vxi + d_r*dxi*momi/mtau;
      vyf = vyi + d_r*dyi*momi/mtau;
      vzf = vzi + d_r*dzi*momi/mtau;
      tf  =  ti + d_r*ei/mtau/(constants::kLightSpeed/(units::centimeter/units::nanosecond));      
    }
  }
  else if (tauproptype=="") {
    idec = 1;
    dxf = dxi; dyf = dyi; dzf = dzi; ef = ei;
    vxf = vxi; vyf = vyi; vzf = vzi; tf = ti;
  }

  LOG("TauPropagation", pDEBUG) << "After decay: " << pdgi << ", E = " << ef << " GeV" << ", idec = " << idec;
  LOG("TauPropagation", pDEBUG) << "  length     = " << sqrt(pow(vxi-vxf,2)+pow(vyi-vyf,2)+pow(vzi-vzf,2)) << " cm";
  LOG("TauPropagation", pDEBUG) << "  Position   = [ " << vxf << " cm, " << vyf << " cm, " << vzf << " cm, " << tf << " ns ]";
  LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";

  vxf *= 1e-2; //from cm(tausic/proposal) to m(geom)
  vyf *= 1e-2;
  vzf *= 1e-2;
  tf  *= 1e-9; //from ns(tausic/proposal) to s(geom)

  double momf = TMath::Sqrt( ef*ef - mtau*mtau );

  std::vector<GHepParticle> PropProd;
  if (idec==1) PropProd = Decay( pdgi, vxf, vyf, vzf, tf, dxf*momf, dyf*momf, dzf*momf, ef );
  else         PropProd.push_back(GHepParticle(pdgi,kIStUndefined,-1,-1,-1,-1,dxf*momf,dyf*momf,dzf*momf,ef,vxf,vyf,vzf,tf));

  return PropProd;

}


std::vector<GHepParticle> TauPropagation::Decay(double pdg, double vx, double vy, double vz, double t, double px, double py, double pz, double e) {

  //decay tau
  TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
  TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( pdg, 1, px, py, pz, e, mtau, -1, -1, -1, -1 );
  Tauola_evt->addParticle(Tauola_tau);
  Tauola::decayOne(Tauola_tau,true,0.,0.,polarization);

  std::vector<GHepParticle> DecProd;
  
  for ( int sec=1; sec<Tauola_evt->getParticleCount(); sec++ ) {
    int spdg   = Tauola_evt->getParticle(sec)->getPdgID();
    double se  = Tauola_evt->getParticle(sec)->getE();
    if ( pdg::IsNeutrino(TMath::Abs(spdg)) ) {
      double spx = Tauola_evt->getParticle(sec)->getPx();
      double spy = Tauola_evt->getParticle(sec)->getPy();
      double spz = Tauola_evt->getParticle(sec)->getPz();
      LOG("TauPropagation", pDEBUG) << "Product: " << spdg << ", E = " << se << " GeV";
      LOG("TauPropagation", pDEBUG) << "  Position   = [ " << vx << " m, " << vy << " m, " << vz << " m, " << t << " s ]";
      LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << spx/se << " , " << spy/se << " , " << spz/se << " ]";
      DecProd.push_back(GHepParticle(spdg,kIStUndefined,-1,-1,-1,-1,spx,spy,spz,se,vx,vy,vz,t));
    }
  }

  delete Tauola_evt;

  return DecProd;

}





