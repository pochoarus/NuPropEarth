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

using namespace std;
using namespace genie;

//const double fREarth_m = 6368e3; //currently we neglect the water layer
const double fREarth_m = 6371e3; //we consider the water/ice layer

namespace genie {
  namespace flux  {

    class IncomingFlux: public GFluxI {

      public :
        IncomingFlux(int pdg, double alpha, double cthmin, double cthmax, double cthmono, double emin, double emax, double emono, double depth, double radius, double height);
        virtual ~IncomingFlux();

        // methods implementing the GENIE GFluxI interface
        virtual PDGCodeList &          FluxParticles (void);
        virtual bool                   GenerateNext  (void);
        virtual double                 MaxEnergy     (void) { return fEmax;      }
        virtual int                    PdgCode       (void) { return fgPdgCI;    }
        virtual double                 Weight        (void) { return 1;          }  
        virtual const TLorentzVector & Momentum      (void) { return fgP4I;      }
        virtual const TLorentzVector & Position      (void) { return fgX4I;      }
        virtual bool                   End           (void) { return false;      }
        virtual long int               Index         (void) { return -1;         }
        virtual void                   Clear            (Option_t * opt)    {}
        virtual void                   GenerateWeighted (bool gen_weighted) {}
        
        // set neutrino to simulate a new interaction (PropMode)
        void           InitNeutrino   (double px, double py, double pz, double e, double x, double y, double z, double t, int pdg, const TLorentzVector & IntVer);
        int            GetPdgI       (void) { return fgPdgCI; }
        TLorentzVector GetP4I        (void) { return fgP4I;   }
        TLorentzVector GetX4I        (void) { return fgX4I;   }

      protected:

        // protected methods
        bool    GenerateNext_1try (void);
        void    Initialize        (void);
        void    ResetSelection    (void);

        // protected data members
        int              fPdg;            
        double           fSpectralIndex;          
        double           fCThmin;          
        double           fCThmax;          
        double           fCThmono;          
        double           fEmin;          
        double           fEmax;          
        double           fEmono;
        double           fDepth;
        double           fRadius;
        double           fHeight;          

        bool             fNewNeutrino;
        TLorentzVector   fgP4I;            
        TLorentzVector   fgX4I;            
        int              fgPdgCI;          
      
        
    };

  } // flux namespace
} // genie namespace

#endif // _IncomingFlux_H_
