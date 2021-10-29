
#include "HadronPropagation.h"

extern "C" { 
  void initialize_tausic_sr_(int *muin);
  void initialize_tausic_sw_(int *muin);
  void tau_transport_sr_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
  void tau_transport_sw_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
}

HadronPropagation::HadronPropagation(ROOTGeomAnalyzer * gd) {

  geom_driver = gd;

  LOG("HadronPropagation", pDEBUG) << "Initializing Random...";
  rnd = RandomGen::Instance();

  fPythia = TPythia6::Instance();

  PdgDB = TDatabasePDG::Instance();

}

void HadronPropagation::ComputeDepth(GHepParticle * p, double &avgrho, double &totlength) { //length(m) depth(g/cm2)

  totlength = 0.;
  avgrho    = 0.;
  std::vector< std::pair<double, const TGeoMaterial*> > MatLengths = geom_driver->ComputeMatLengths(*p->X4(),*p->P4());
  for ( auto sitr : MatLengths ) {
    double length  = sitr.first; 
    double rho     = sitr.second->GetDensity();
    totlength += length;
    avgrho    += length*rho;
  }
  avgrho = avgrho/totlength;
  LOG("HadronPropagation", pDEBUG) << "Sum ---> Avg. Rho = " << avgrho << " g/cm^3 ; Length = " << totlength << " m";

}


std::vector<GHepParticle> HadronPropagation::Propagate(GHepParticle * hadron, string pclass) { 

  int pdg         = hadron->Pdg();
  double kf       = fPythia->Pycomp(pdg);
  double lifetime = fPythia->GetPMAS(kf,4)*1e-3; //pythia gives lifetime in mm (c*t)
  double mass     = fPythia->GetPMAS(kf,1); //GeV

  double vx = hadron->Vx();
  double vy = hadron->Vy();
  double vz = hadron->Vz();
  double t  = hadron->Vt();

  double mom = hadron->P4()->P();
  double dx  = hadron->Px()/mom;
  double dy  = hadron->Py()/mom;
  double dz  = hadron->Pz()/mom;
  double e   = TMath::Sqrt( mom*mom + mass*mass );
  
  LOG("HadronPropagation", pDEBUG) << "Before decay: " << pdg << ", E = " << e << " GeV";
  LOG("HadronPropagation", pDEBUG) << "  Position   = [ " << vx << " m, " << vy << " m, " << vz << " m, " << t << " s ]";
  LOG("HadronPropagation", pDEBUG) << "  Direction  = [ " << dx << " , " << dy << " , " << dz << " ]";

  int idec = 0; //decay flag (0=not decay // 1=decay // 2=destroyed)

  double avgrho,depthi; //g/cm3, m
  ComputeDepth(hadron,avgrho,depthi);

  double depth0 = depthi;

  while( depthi>1 ) {

    // length before it decays [m]
    double lengthdec = -TMath::Log(rnd->RndGen().Rndm()) * mom * lifetime / mass;

    // length before it interact [m]
    double lengthint = -TMath::Log(rnd->RndEvg().Rndm()) * HadronInteraction(pclass,mass,e,avgrho);

    LOG("HadronPropagation", pDEBUG) << e << "  " << HadronInteraction(pclass,mass,e,avgrho) << "  " << mom * lifetime / mass;

    //particle to close to the surface
    if (lengthdec>depthi && lengthint>depthi) {
      vx += depthi*dx;
      vy += depthi*dy;
      vz += depthi*dz;
      t  += depthi/(constants::kLightSpeed/(units::meter/units::second));      
      break;
    }

    if (lengthint>lengthdec) {
      idec = 1; //particle decayed 
      vx += lengthdec*dx;
      vy += lengthdec*dy;
      vz += lengthdec*dz;
      t  += lengthdec/(constants::kLightSpeed/(units::meter/units::second));      
      break;
    }

    double z = HadronInelasticity(pclass);
    LOG("HadronPropagation", pDEBUG) << "Interaction -> lenght =  " << lengthint << ", z = " << z;

    e *= z;
    if ( e<mass ) {
      idec = 2; //particle dissapear      
      break; 
    }

    vx += lengthint*dx;
    vy += lengthint*dy;
    vz += lengthint*dz;
    t  += lengthint/(constants::kLightSpeed/(units::meter/units::second));      

    mom = TMath::Sqrt( e*e - mass*mass );
    hadron->SetMomentum(dx*mom, dy*mom, dz*mom, e);
    hadron->SetPosition(vx,vy,vz,t);
    ComputeDepth(hadron,avgrho,depthi);

  }

  LOG("HadronPropagation", pDEBUG) << "After decay: " << pdg << ", E = " << e << " GeV" << ", idec = " << idec;
  LOG("HadronPropagation", pDEBUG) << "depth0 = " << depth0 << "  depthi = " << depthi << " -> " << depth0-depthi; 
  LOG("HadronPropagation", pDEBUG) << "  Position   = [ " << vx << " m, " << vy << " m, " << vz << " m, " << t << " s ]";
  LOG("HadronPropagation", pDEBUG) << "  Direction  = [ " << dx << " , " << dy << " , " << dz << " ]";


  std::vector<GHepParticle> PropProd;
  if      (idec==1) PropProd = Decay( pdg, kf, vx, vy, vz, t, dx, dy, dz, e );
  else if (idec==0) {
    LOG("HadronPropagation", pDEBUG) << "Reached the surface";  
    PropProd.push_back(GHepParticle(pdg,kIStUndefined,-1,-1,-1,-1,dx*mom,dy*mom,dz*mom,e,vx,vy,vz,t));
  }

  return PropProd;

}


std::vector<GHepParticle> HadronPropagation::Decay(double pdg, double kf, double vx, double vy, double vz, double t, double dx, double dy, double dz, double e) {

  fPythia->Py1ent( -1, pdg, e, acos(dz), atan2(dy,dx) );

  int decflag = fPythia->GetMDCY(kf,1);
  fPythia->SetMDCY(kf,1,1); //enable the decay of this particular hadron          
  fPythia->Pyexec();
  fPythia->SetMDCY(kf,1,decflag); //go back to original

  fPythia->GetPrimaries();
  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");

  std::vector<GHepParticle> DecProd;

  TMCParticle * particle = 0;
  TIter piter(pythia_particles);      
  while( (particle = (TMCParticle *) piter.Next()) ) {
    if (particle->GetKS()!=1) continue;
    int spdg   = particle->GetKF();
    double se  = particle->GetEnergy();
    
    bool fill = false;
    if ( pdg::IsNeutrino(TMath::Abs(spdg)) || pdg::IsTau(TMath::Abs(spdg)) ) fill = true;
    else {
      TParticlePDG * partinfo = PdgDB->GetParticle(spdg);
      string pclass = partinfo->ParticleClass();
      if( pclass=="CharmedMeson" || pclass=="CharmedBaryon" ) fill = true;
    }

    if (fill) {
      double spx = particle->GetPx();
      double spy = particle->GetPy();
      double spz = particle->GetPz();
      LOG("HadronPropagation", pDEBUG) << "Product: " << spdg << ", E = " << se << " GeV";
      LOG("HadronPropagation", pDEBUG) << "  Position   = [ " << vx << " m, " << vy << " m, " << vz << " m, " << t << " s ]";
      LOG("HadronPropagation", pDEBUG) << "  Direction  = [ " << spx/se << " , " << spy/se << " , " << spz/se << " ]";
      DecProd.push_back(GHepParticle(spdg,kIStUndefined,-1,-1,-1,-1,spx,spy,spz,se,vx,vy,vz,t));
    }
  }
  delete particle;
  pythia_particles->Clear("C");

  return DecProd;

}

//___________________________________________________________________________
double HadronPropagation::HadronInteraction(string pclass, double mass, double E, double rho) {

  double xsec_p = 0.; //mbarn

  double p = TMath::Sqrt(E*E-mass*mass);

  // from https://inspirehep.net/literature/1707889 (section 5.2.3)
  // from https://inspirehep.net/literature/469835 (section C.1)
  if (pclass=="CharmedMeson") {
    if      (E<1e1) xsec_p = 0.;
    else if (E<1e3) xsec_p = 0.87549*(12.3-7.7*TMath::Power(p,-2.12)+0.0326*TMath::Power(TMath::Log(p),2)+0.738*TMath::Log(p));
    else if (E<1e6) xsec_p = 0.87549*(11.01*TMath::Power(p*p,0.0395));
    else            xsec_p = TMath::Exp(1.891+0.2095*TMath::Log10(E))+1.263*TMath::Log10(E)-2.157;
  }
  else if (pclass=="CharmedBaryon") {
    if      (E<1e1) xsec_p = 0.;
    if      (E<1e3) xsec_p = 1.22442*(12.3-7.7*TMath::Power(p,-2.12)+0.0326*TMath::Power(TMath::Log(p),2)+0.738*TMath::Log(p));
    else if (E<1e6) xsec_p = 1.22442*(11.01*TMath::Power(p*p,0.0395));
    else            xsec_p = TMath::Exp(2.269+0.207*TMath::Log10(E))+1.277*TMath::Log10(E)-0.9907;
  }
  else if (pclass=="B-Meson") {
    if      (E<1e1) xsec_p = 0.;
    if      (E<1e3) xsec_p = 0.78335*(12.3-7.7*TMath::Power(p,-2.12)+0.0326*TMath::Power(TMath::Log(p),2)+0.738*TMath::Log(p));
    else if (E<1e6) xsec_p = 0.78335*(11.01*TMath::Power(p*p,0.0395));
    else            xsec_p = TMath::Exp(1.851+0.2094*TMath::Log10(E))+0.7279*TMath::Log10(E)-1.042;
  }
  else if (pclass=="B-Baryon") {
    if      (E<1e1) xsec_p = 0.;
    if      (E<1e3) xsec_p = 1.15313*(12.3-7.7*TMath::Power(p,-2.12)+0.0326*TMath::Power(TMath::Log(p),2)+0.738*TMath::Log(p));
    else if (E<1e6) xsec_p = 1.15313*(11.01*TMath::Power(p*p,0.0395));
    else            xsec_p = TMath::Exp(2.230+0.207*TMath::Log10(E))+1.086*TMath::Log10(E)-0.9026;
  }

  double s45 = (xsec_p-45.)/30;  
  double xsec_o = (1.-4.*TMath::Power(s45,2))*342.6 + s45*(2.*s45-1.)*271.6 + s45*(2.*s45+1.)*399.7; //mbarn

  double xsec_h2o = (0.11*xsec_p/1.00784 + 0.89*xsec_o/15.999); //mbarn

  return 1e25 / TMath::Na() / rho / xsec_h2o; //m 

}
//___________________________________________________________________________
double HadronPropagation::HadronInelasticity(string pclass) {

  double y = 0.;
  
  // from https://inspirehep.net/literature/1707889 (section 5.2.3)
  if      (pclass=="CharmedMeson" ) y = rnd->RndGen().Gaus(0.56,0.2);
  else if (pclass=="CharmedBaryon") y = rnd->RndGen().Gaus(0.59,0.2);
  else if (pclass=="B-Meson")       y = rnd->RndGen().Gaus(0.80,0.2);
  else if (pclass=="B-Baryon")      y = rnd->RndGen().Gaus(0.70,0.2);
  
  if      (y<0) y=0;
  else if (y>1) y=1;

  return y;

}





