
#include "TauPropagation.h"

extern "C" { 
  void initialize_tausic_sr_(int *muin);
  void tau_transport_sr_(double *VX,double *VY,double *VZ,double *CX,double *CY,double *CZ,double *E,double *Depth,double *T,double *Rho,int *flag1,int *flag2,double *VXF,double *VYF,double *VZF,double *CXF,double *CYF,double *CZF,double *EF,double *DepthF,double *TF, int *IDEC, int *ITFLAG, double *TAUTI, double *TAUTF);
}

TauPropagation::TauPropagation(string ptype, int seed) {

  tauproptype = ptype;

  LOG("TauPropagation", pDEBUG) << "Initializing TAUOLA...";
  Tauola::initialize();
  Tauola::setSeed(seed,seed,seed);
  mtau         = Tauola::getTauMass();        //use this mass to avoid energy conservation warning in tauola
  d_lifetime   = Tauola::tau_lifetime * 1e-1; //lifetime from mm(tauola) to cm
  polarization = 1;                           //tau-(P=-1) & tau+(P=-1) however in tauola its is fliped

  if (tauproptype=="TAUSIC") {
    TAUIN = 2;
    TAUMODEL = 1;
    ITFLAG = 0; //set lifetime inside tausic
    LOG("TauPropagation", pDEBUG) << "Initializing TAUSIC...";
    string tausicpath = gSystem->Getenv("TAUSIC");
    string cp_command = "cp "+tausicpath+"/tausic*.dat .";
    system(cp_command.c_str());
    initialize_tausic_sr_(&TAUIN);      
    system("rm ./tausic*.dat");
  }

  LOG("TauPropagation", pDEBUG) << "Initializing Random...";
  rnd = RandomGen::Instance();

}


vector<GHepParticle> TauPropagation::Propagate(GHepParticle * tau, double lengthi, double depthi) { //length(cm) depth(g/cm2)

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

  if      (tauproptype=="TAUSIC") {
    depthi = depthi/rockdensity;
    double tauti,tautf; //dummy only used when ITFLAG!=0
    double depthf; //total distance travel after propagation
    tau_transport_sr_(&vxi,&vyi,&vzi,&dxi,&dyi,&dzi,&ei,&depthi,&ti,&rockdensity,&TAUMODEL,&TAUMODEL,&vxf,&vyf,&vzf,&dxf,&dyf,&dzf,&ef,&depthf,&tf,&idec,&ITFLAG,&tauti,&tautf);
  }
  else if (tauproptype=="NOELOSS") {
    dxf = dxi; dyf = dyi; dzf = dzi; ef = ei;
    // position based on lifetime
    double d_r = -TMath::Log( rnd->RndGen().Rndm() ) * d_lifetime;    
    if ( d_r*momi/mtau>lengthi ) {
      idec = 0;
      vxf = vxi + lengthi*dxi;
      vyf = vyi + lengthi*dyi;
      vzf = vzi + lengthi*dzi;
      tf  =  ti + lengthi/(constants::kLightSpeed/(units::centimeter/units::nanosecond));      
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

  LOG("TauPropagation", pDEBUG) << "After decay: " << pdgi << ", E = " << ef << " GeV";
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



