
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
  Tauola::setSeed(seed,seed,seed);
  mtau         = Tauola::getTauMass();        //use this mass to avoid energy conservation warning in tauola
  d_lifetime   = Tauola::tau_lifetime/(constants::kLightSpeed/(units::millimeter/units::nanosecond)); //lifetime in ns
  polarization = 1;                           //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped
  LOG("TauPropagation", pDEBUG) << d_lifetime;


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

  LOG("TauPropagation", pDEBUG) << "Initializing Random...";
  rnd = RandomGen::Instance();

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


vector<GHepParticle> TauPropagation::Propagate(GHepParticle * tau) { 

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

  int idec; //decay flag (0=not decay // >0=decay)
  double vxf,vyf,vzf,tf; //position after propagation in cm/ns
  double dxf,dyf,dzf,ef; //direction after propagation

  double avgrho,depthi;
  ComputeDepth(tau,avgrho,depthi);

  double tauti = -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime; //ns
  LOG("TauPropagation", pDEBUG) << tauti << " ns";

  if      (tauproptype=="TAUSIC") {
    while ( depthi>1 ) { //threshold in 1cm
      double tautf; //remaining lifetime [ns]
      double depthf; //total distance travel after propagation [cm]
      if   ( avgrho>2. ) tau_transport_sr_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&avgrho,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);  
      else               tau_transport_sw_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&avgrho,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);
      LOG("TauPropagation", pDEBUG) << "  depthf  = " << depthf << " cm";

      if ( idec==1 ) break; //tau has decayed

      double momf = TMath::Sqrt( ef*ef - mtau*mtau );
      tau->SetMomentum(dxf*momf, dyf*momf, dzf*momf, ef);
      tau->SetPosition(vxf*1e-2,vyf*1e-2,vzf*1e-2,tf*1e-9); //from cm [tausic] to m [geom]
      ComputeDepth(tau,avgrho,depthi);
      
      tauti = tautf;
      vxi = vxf; vyi = vyf; vzi = vzf; ti = tf;
      dxi = dxf; dyi = dyf; dzi = dzf; ei = ef; 
      
      LOG("TauPropagation", pDEBUG) << tauti << " ns";
    }

  }
  else if (tauproptype=="NOELOSS") {
    dxf = dxi; dyf = dyi; dzf = dzi; ef = ei;
    double d_r  = tauti*(constants::kLightSpeed/(units::centimeter/units::nanosecond)); //cm
    LOG("TauPropagation", pDEBUG) << d_r;
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
  LOG("TauPropagation", pDEBUG) << "  Position   = [ " << vxf << " cm, " << vyf << " cm, " << vzf << " cm, " << tf << " ns ]";
  LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << dxf << " , " << dyf << " , " << dzf << " ]";

  vxf *= 1e-2; //from cm(tausic) to m(geom)
  vyf *= 1e-2;
  vzf *= 1e-2;
  tf  *= 1e-9; //from ns(tausic) to s(geom)

  double momf = TMath::Sqrt( ef*ef - mtau*mtau );

  vector<GHepParticle> PropProd;
  if (idec>0) PropProd = Decay( pdgi, vxf, vyf, vzf, tf, dxf*momf, dyf*momf, dzf*momf, ef );
  else        PropProd.push_back(GHepParticle(pdgi,kIStUndefined,-1,-1,-1,-1,dxf*momf,dyf*momf,dzf*momf,ef,vxf,vyf,vzf,tf));

  return PropProd;

}


vector<GHepParticle> TauPropagation::Decay(double pdg, double vx, double vy, double vz, double t, double px, double py, double pz, double e) {

  //decay tau
  TauolaHEPEVTEvent * Tauola_evt = new TauolaHEPEVTEvent();
  TauolaHEPEVTParticle *Tauola_tau = new TauolaHEPEVTParticle( pdg, 1, px, py, pz, e, mtau, -1, -1, -1, -1 );
  Tauola_evt->addParticle(Tauola_tau);
  Tauola::decayOne(Tauola_tau,true,0.,0.,polarization);

  vector<GHepParticle> DecProd;
  
  for ( int sec=1; sec<Tauola_evt->getParticleCount(); sec++ ) {
    int spdg   = Tauola_evt->getParticle(sec)->getPdgID();
    double se  = Tauola_evt->getParticle(sec)->getE();
    if ( pdg::IsNeutrino(TMath::Abs(spdg)) ) {
      double spx = Tauola_evt->getParticle(sec)->getPx();
      double spy = Tauola_evt->getParticle(sec)->getPy();
      double spz = Tauola_evt->getParticle(sec)->getPz();
      LOG("TauPropagation", pDEBUG) << "Product: " << spdg << ", E = " << se << " GeV";
      LOG("TauPropagation", pDEBUG) << "  Position   = [ " << vx << " , m" << vy << " , m" << vz << " , m" << t << " s ]";
      LOG("TauPropagation", pDEBUG) << "  Direction  = [ " << spx/se << " , " << spy/se << " , " << spz/se << " ]";
      DecProd.push_back(GHepParticle(spdg,kIStUndefined,-1,-1,-1,-1,spx,spy,spz,se,vx,vy,vz,t));
    }
  }

  delete Tauola_evt;

  return DecProd;

}



