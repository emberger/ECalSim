
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

#include "Minuit2/MnUserParameters.h"
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

void findShowercenter(Int_t minevent, Int_t maxevent);
Double_t FindClosestApproach();
void CalculateDeviation();

void CalcCOGPion(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent, Int_t photonNR=1);

void CalcCOG(Int_t minevent, Int_t maxevent);
/////////////////////////////////////////////////////////////////////////////////////////////
void FitCOGs( Int_t minevent, Int_t maxevent, Double_t tileLen);

void CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);
//////////////////////////////////////////////////////////////////////////////////////////////

void CalcCOGwithFit(Int_t minevent, Int_t maxevent);

void FitCOGsPion( Int_t minevent, Int_t maxevent, Double_t tileLen, Bool_t isPion=false, Int_t photonNR=1);

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
Double_t Eges;        //Energy in event
Int_t showerstart;

Double_t histsizeX=100;
Double_t histsizeY=100;
Double_t histsizeZ=50;

std::string savepath;

Bool_t pathset=false;

std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > coglist;
std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > COGCollection;
std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParamsGamma;
std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > FitParams;//to be deleted, needed for old analysis
std::vector<std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t> > > FitParamsPions;
std::vector<std::tuple<Double_t, Double_t> > ClosestApproach;

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

TH3D * h = new TH3D("ECalEvent","ECalEvent",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);

TH3D * hA = new TH3D("ECalEvent1","ECalEvent1",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
TH3D * hB = new TH3D("ECalEvent2","ECalEvent2",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);                                //plot Event

TH2D * h1 = new TH2D("h1", "h1", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);
//  clustering
TH2D * h2 = new TH2D("h2", "h2", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);                                 //  clustering

TH3D * h3 = new TH3D("h3", "h3",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);           //  plot cogs

TH1D * delx = new TH1D("DeltaX", "DeltaX", 1000,-500,500);
TH1D * dely = new TH1D("DeltaY", "DeltaY", 1000,-500,500);
TH1D * delz = new TH1D("DeltaZ", "DeltaZ", 800,-500,1500);




};
#endif
