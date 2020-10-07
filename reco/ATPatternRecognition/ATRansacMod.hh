/*******************************************************************
// Basic RANSAC Class                                              *
// Author: J.C. Zamora, jczamorac@gmail.com                        *
// University of Sao Paulo, 26-08-2020                             *
********************************************************************/

#ifndef ATRANSACMOD_H
#define ATRANSACMOD_H

#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>

// FairRoot classes
#include "FairRootManager.h"
#include "FairLogger.h"

// ROOT Headers
#include <TObject.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>
#include "TApplication.h"
#include "TROOT.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ATHit.hh"
#include "ATEvent.hh"
#include "ATProtoEvent.hh"
#include "ATProtoQuadrant.hh"
#include "ATDigiPar.hh"
#include "AtTpcMap.h"
#include "ATTrack.hh"

using namespace std;


#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"



class ATRansacMod : public TObject
{

  public:
      ATRansacMod();
      ~ATRansacMod();


      void Reset();
	    void Init(ATEvent *event);
	    void Solve();
      vector<int> RandSam(vector<int> indX, Int_t mode);
      void EstimModel(const std::vector<int>  samplesIdx);
      double EstimError(int i);
      void CalcRANSACMod(ATEvent *event);
	    vector<double> GetChargeOfTracks();
	    vector<double> GetTrackLength();
      vector<double> GetPDF(const std::vector<int>  samplesIdx);
      void SetAvCharge(double charge){Avcharge = charge;};
      void SetRanSamMode(Int_t mode){fRandSamplMode = mode;};
      double GetAvCharge(){return Avcharge;};
	    double Fit3D(vector<int> inliners, TVector3& V1, TVector3& V2);
      void SetDistanceThreshold(Float_t threshold) { fRANSACThreshold = threshold;};
      void SetMinHitsLine(Int_t nhits) { fRANSACMinPoints = nhits;};
      void SetNumItera(Int_t niterations) { fRANSACMaxIteration = niterations;};
      TVector3 GetVertex1(){return fVertex_1;};
      TVector3 GetVertex2(){return fVertex_2;};
      Double_t GetVertexTime(){return fVertexTime;};
      TVector3 GetVertexMean(){return fVertex_mean;};
      std::vector<ATTrack> GetTrackCand(){return fTrackCand;};
      TVector3 ClosestPoint2Lines(TVector3 d1, TVector3 pt1, TVector3 d2, TVector3 pt2);
      void SetChargeThres(double value){fChargeThres = value;};

      struct Cluster // return type of structure
        {
          double ClusterStrength;		// strength
          size_t ClusterSize;			// size
          double ClusterChi2;			// Chi2
          std::vector<int> ClusterIndex;			// Indices
          TVector3 ClusterFitP1;			// point 1 from the fitted line
          TVector3 ClusterFitP2;			// point 2 from the fitted line
        };


      typedef std::vector<Cluster> AllClusters;
      std::vector<ATTrack*> Clusters2Tracks( AllClusters NClusters, ATEvent *event);

      void SetCluster(const std::vector<int> samplesIdx, const double cost, const double Chi2, TVector3 CP1, TVector3 CP2);
      inline AllClusters GetClusters(){return cluster_vector;}
      AllClusters cluster_vector;




  protected:

      void FindVertex(std::vector<ATTrack*> tracks);
      void FindVertexOneTrack(std::vector<ATTrack*> tracks);


      TVector3 fVertex_1;
      TVector3 fVertex_2;
      TVector3 fVertex_mean;
      Double_t fVertexTime;
      Double_t fMinimum;
      std::vector<ATTrack> fTrackCand; //Candidate tracks
      std::pair<Int_t,Int_t> fVertex_tracks; // ID of the tracks that form the best vertex
      Int_t fLineDistThreshold;
      Int_t fRandSamplMode;

      vector<double> vX, vY, vZ, vQ;
	    vector<double> vTrackCharge;
      float fRANSACMinPoints;
	    float fRANSACPointThreshold;
	    float fRANSACChargeThreshold;
	    float fRANSACThreshold;
	    float fRANSACMaxIteration;
	    int fNumberOfTracksMax;
	    int fOriginalCloudSize;
	    double fTotalCharge;
	    int fVerbose;
      double Avcharge;
      double fChargeThres;



  public:

    TVector3 Vs;
    TVector3 Ps;


  ClassDef(ATRansacMod, 1);


};

#endif
