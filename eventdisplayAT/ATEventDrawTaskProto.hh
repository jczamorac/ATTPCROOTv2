#ifndef ATEVENTDRAWTASKPROTO_H
#define ATEVENTDRAWTASKPROTO_H


// FairRoot classes
#include "FairTask.h"
#include "FairLogger.h"

// ROOT classes
#include "TEvePointSet.h"
#include "TEveGeoShape.h"
#include "TEveBoxSet.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TPaletteAxis.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TH3.h"
#include "TGraph.h"
#include "TH2Poly.h"

#include "ATEventManagerProto.hh"
//#include "ATRiemannTrack.hh"
//#include "ATRiemannHit.hh"
#include "ATRawEvent.hh"
#include "ATEvent.hh"
#include "ATProtoEvent.hh"
#include "ATHoughSpaceLine.hh"
#include "ATHoughSpaceCircle.hh"
#include "ATHoughSpace.hh"
#include "ATHit.hh"
#include "AtTpcMap.h"
#include "ATProtoQuadrant.hh"
#include <fstream>

#ifndef __CINT__ // Boost
#include <boost/multi_array.hpp>
#endif //__CINT__

class ATEventDrawTaskProto : public FairTask
{
  public :
    ATEventDrawTaskProto();

    virtual ~ATEventDrawTaskProto();

    virtual InitStatus Init();
    virtual void Exec(Option_t* option);
    void Reset();
    void SetHitAttributes(Color_t, Size_t, Style_t);
    void Set3DHitStyleBar();
    void Set3DHitStyleBox();

    static void SelectPad(const char *rawevt);

    void SetProtoMap(TString map) {fMap = map;}

  protected:
    virtual void DrawPadWave();
    virtual void DrawPadPlane();
    virtual void DrawPadAll();
    virtual void DrawMesh();
    virtual void DrawProtoSpace();
    virtual void DrawProtoEL();
    virtual void DrawProtoHough();

    void DrawHitPoints();
    void DrawProtoPattern();

    void UpdateCvsPadWave();
    void UpdateCvsPadPlane();
    void UpdateCvsPadAll();
    void UpdateCvsMesh();
    void UpdateCvsProtoQ();
    void UpdateCvsProtoEL();


     //Basic types

    Int_t fMultiHit;
    Int_t f3DHitStyle;
    Bool_t fUnpackHough;
    Bool_t fIsCircularHough;
    Bool_t fIsLinearHough;
    Bool_t fIsRawData;
    Color_t fHitColor;
    Size_t  fHitSize;
    Style_t fHitStyle;

    // ROOT Objects
    TPaletteAxis *fPadPlanePal;


    TH3F* f3DHist;

    TH1I*    fPadAll[300];
    TH1D*    fPhiDistr[5];
    TH1I*    fPadWave;
    TH2Poly* fPadPlane;
    TH1F*    fMesh;
    TH2F*    fQuadrant1;
    TH2F*    fQuadrant2;
    TH2F*    fQuadrant3;
    TH2F*    fQuadrant4;
    TGraph*  fQHitPattern[4];
    TGraph*  fQELossPattern[4];
    TF1*     fHoughFit[4];

    TCanvas* fCvsPadWave;
    TCanvas* fCvsPadPlane;
    TCanvas* fCvsPadAll;
    TCanvas* fCvsMesh;
    TCanvas* fCvsQuadrant1;
    TCanvas* fCvsQuadrant2;
    TCanvas* fCvsQuadrant3;
    TCanvas* fCvsQuadrant4;
    TCanvas* fCvsELQuadrant1;
    TCanvas* fCvsELQuadrant2;
    TCanvas* fCvsELQuadrant3;
    TCanvas* fCvsELQuadrant4;


    TF1 *fHoughLinearFit;
    TString fMap;

    TClonesArray* fHitArray;
    TClonesArray* fRawEventArray;
    TClonesArray* fHoughSpaceArray;
    TClonesArray* fProtoEventArray;

    TEvePointSet* fHitSet;
    TEveBoxSet* fhitBoxSet;

    /// ATTPCROOT objects

    ATEventManagerProto*  fEventManager;
    ATHit*                fIniHit;
    AtTpcMap*             fDetmap;
    ATRawEvent*           fRawevent;
    ATHoughSpaceLine*     fHoughSpaceLine;




    ClassDef(ATEventDrawTaskProto,1);
  };

  #endif