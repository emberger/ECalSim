
#ifndef TROOTAnalysis_hh
#define TROOTAnalysis_hh


#include "B4ROOTEvent.hh"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TBrowser.h"
#include "TColor.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraph2D.h"
#include "TAttMarker.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TImage.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minimizer.hh"
#include "Minuit2/Minuit2Minimizer.h"

#include <chrono>
#include <string>
#include <cstring>

using namespace ROOT::Math;
using namespace ROOT::Minuit2;



class TROOTAnalysis {

public:
TROOTAnalysis(TChain* ch);
~TROOTAnalysis();


void PrintERes();
void PlotRMSx();
void plotEvent(Int_t pev);
void plotEventPion(Int_t pev);
void plotCOGs();
void PlotProjection(Int_t distance);
void PlotChimap(Int_t event);


Double_t FindClosestApproach();
void CalculateDeviation();
void GetInvariantMass();
//void MinimizeClosestApproach();
void PionLocator();

void CalcCOGPion(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent, Int_t photonNR=1);

void CalcCOG(Int_t minevent, Int_t maxevent);
/////////////////////////////////////////////////////////////////////////////////////////////
void FitCOGs( Int_t minevent, Int_t maxevent);

void CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);
//////////////////////////////////////////////////////////////////////////////////////////////



void FitCOGsPion( Int_t minevent, Int_t maxevent, Bool_t isPion=false, Int_t photonNR=1);

void AnalyzePions(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);
//  void CalcCOGwithFit(Int_t minlayer, Int_t maxlayer);
void CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);


void PrintFitHists(Int_t minevent, Int_t maxevent);



void SetPath(std::string path){
        savepath=path;
        pathset=true;
}


private:
Int_t nofEntries;     // number of events in Tree

Double_t Eges;  //Energy in event

Double_t EPhot1;
Double_t EPhot2;



Double_t tiledimX;
Double_t tiledimY;
Double_t calsizeXY;
Double_t AbsoThickness;
Double_t GapThickness;
Double_t nofLayers;

Double_t histsizeX;
Double_t histsizeY;
Double_t histsizeZ;

std::string savepath;

Bool_t pathset=false;

std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > coglist;
std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > COGCollection;

std::vector<std::tuple<Double_t, Double_t, Double_t> > showerCOGPhoton1;
std::vector<std::tuple<Double_t, Double_t, Double_t> > showerCOGPhoton2;


std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParamsGamma;
std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParams;//to be deleted, needed for old analysis
std::vector<std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > > FitParamsPions;

std::vector<Double_t > EnergyPhoton1;
std::vector<Double_t > EnergyPhoton2;
std::vector<Double_t> InvariantMass;

std::vector<std::tuple<Double_t,Double_t,Double_t, Double_t> > ClosestApproach;

std::vector<std::tuple<Double_t, Double_t, Double_t> > DeviationFromGun;

std::vector<std::tuple<Double_t, Double_t, Double_t> > ClusteredHits;
std::vector<Double_t> showerCenters;

TTree* EcalTree;
B4ROOTEvent * Cevent;

TCanvas * c2 = new TCanvas("COGs", "COGs");

TH1D * fitX = new TH1D("X Fit","X Fit", 4000,-600,600);
TH1D * fitY = new TH1D("Y Fit","Y Fit", 4000,-600,600);

TH1D * fitSX = new TH1D("X Slope Fit","X Slope Fit", 800,-2,2);
TH1D * fitSY = new TH1D("Y Slope Fit","Y Slope Fit", 800,-2,2);


TH1D * er1= new TH1D("errx", "errx", 50,0,50);
TH1D * er2= new TH1D("erry", "erry", 50,0,50);

TH3D * h;

TH3D * hA;
TH3D * hB;

TH2D * h1;
//  clustering
TH2D * h2;

TH3D * h3;

TH1D * delx = new TH1D("DeltaX", "DeltaX", 1000,-500,500);
TH1D * dely = new TH1D("DeltaY", "DeltaY", 1000,-500,500);
TH1D * delz = new TH1D("DeltaZ", "DeltaZ", 800,-1000,1000);
TH1D * appdist1=new TH1D("Closest Approach", "ClosestApproach", 500,0,500);




};
#endif
