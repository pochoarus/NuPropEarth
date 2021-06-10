#ifndef _IncomingFlux_H_
#define _IncomingFlux_H_

#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

#include <TLorentzVector.h>
#include <TRotation.h>

#include "Framework/EventGen/GFluxI.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/GHEP/GHepParticle.h"

using namespace std;
using namespace genie;

const double fREarth_m = 6371e3; //we consider the water/ice layer

namespace genie {
  namespace flux  {

    class IncomingFlux: public GFluxI {

      public :
        IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double emin, double emax, double detpos[3], double radius, double height);
        virtual ~IncomingFlux();

        // methods implementing the GENIE GFluxI interface
        virtual PDGCodeList &          FluxParticles (void);
        virtual bool                   GenerateNext  (void);
        virtual double                 MaxEnergy     (void) { return fEmax;            }
        virtual int                    PdgCode       (void) { return fNeutrino->Pdg(); }
        virtual double                 Weight        (void) { return 1;                }  
        virtual const TLorentzVector & Momentum      (void) { return *fNeutrino->P4(); }
        virtual const TLorentzVector & Position      (void) { return *fNeutrino->X4(); }
        virtual bool                   End           (void) { return false;            }
        virtual long int               Index         (void) { return -1;               }
        virtual void                   Clear            (Option_t * opt)    {}
        virtual void                   GenerateWeighted (bool gen_weighted) {}
        
        // set neutrino to simulate a new interaction (PropMode)
        void           InitNeutrino (GHepParticle Nu);
        GHepParticle * GetNeutrino  (void) { return fNeutrino; }

      protected:

        // protected methods
        bool    GenerateNext_1try (void);

        // protected data members
        PDGCodeList      fPdg;
        double           fSpectralIndex;          
        double           fCThmin;          
        double           fCThmax;          
        double           fEmin;          
        double           fEmax;          
        double           fDetPos[3];
        double           fRadius;
        double           fHeight;          

        bool             fNewNeutrino;
        GHepParticle *   fNeutrino;      
      
        
    };

  } // flux namespace
} // genie namespace

#endif // _IncomingFlux_H_
