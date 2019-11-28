#ifndef _GENIE_TR_JOB_DRIVER_H_
#define _GENIE_TR_JOB_DRIVER_H_

#include <string>
#include <map>

#include <TLorentzVector.h>
#include <TBits.h>
#include <TGeoMaterial.h>

#include "Framework/ParticleData/PDGCodeList.h"

using std::string;
using std::map;

namespace genie {

class EventRecord;
class GFluxI;
class GeomAnalyzerI;
class GENIE;
class GEVGPool;

class GTRJDriver {

public :
  GTRJDriver();
 ~GTRJDriver();

  // configure TR job
  void SetEventGeneratorList       (string listname);
  void UseFluxDriver               (GFluxI * flux);
  void UseGeomAnalyzer             (GeomAnalyzerI * geom);
  void Configure                   ();

  // generate single neutrino event for input flux & geometry
  EventRecord * GenerateEvent (void);

private:
 
  // private methods:
  bool          ComputeInteraction              (void);
  void          GenerateEventKinematics         (void);
  double        EvalXsec                        (int mpdg, int nupdg, double Enu);

  // private data members:
  GEVGPool *      fGPool;              ///< A pool of GEVGDrivers properly configured event generation drivers / one per init state
  GFluxI *        fFluxDriver;         ///< [input] neutrino flux driver
  GeomAnalyzerI * fGeomAnalyzer;       ///< [input] detector geometry analyzer
  double          fCurdL;              ///< [current] distance betwen origin and interaction vertex
  EventRecord *   fCurEvt;             ///< [current] generated event
  int             fCurTgtPdg;          ///< [current] selected target material PDG code
  string          fEventGenList;       ///< [config] list of event generators loaded by this driver (what used to be the $GEVGL setting)
  TBits *         fUnphysEventMask;    ///< [config] controls whether unphysical events are returned (what used to be the $GUNPHYSMASK setting)

};

}      // genie namespace
#endif // _GENIE_TR_JOB_DRIVER_H_