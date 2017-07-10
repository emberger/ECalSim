#include "B4ROOTEvent.hh"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TAttMarker.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minimizer.hh"
#include "Minuit2/Minuit2Minimizer.h"



#include <chrono>
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

class TROOTAnalysis{

public:
  TROOTAnalysis(TChain* ch);
  ~TROOTAnalysis();

  void plotEvent(Int_t pev);
  void plotCOGs();

  void CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);
  void CalcCOG(Int_t minevent, Int_t maxevent);

  void CalcCOGwithFit(Int_t minevent, Int_t maxevent);

  void FitCOGs( Int_t minevent, Int_t maxevent);

//  void CalcCOGwithFit(Int_t minlayer, Int_t maxlayer);
  void CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent);


  void PrintFitHists();


private:
  Int_t nofEntries;   // number of events in Tree
  Double_t Eges;      //Energy in event
  Int_t showerstart;

  TTree* EcalTree;
  B4ROOTEvent * Cevent;

  std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t>> coglist;
  std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t>>> COGCollection;
  std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t>> FitParams;
  //std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t, Double_t, Double_t>> FitParams;

  std::vector<std::tuple<Double_t, Double_t, Double_t>> ClusteredHits;





  TH1D * er1= new TH1D("errx", "errx", 50,0,50);
  TH1D * er2= new TH1D("erry", "erry", 50,0,50);

  TCanvas * c2 = new TCanvas("COGs", "COGs");

  Double_t histsizeX=100;
  Double_t histsizeY=100;
  Double_t histsizeZ=50;

  TH2D * h1 = new TH2D("h1", "h1", histsizeX,0,histsizeX,histsizeY,0,histsizeY);
  TH2D * h2 = new TH2D("h2", "h2", histsizeX,0,histsizeX,histsizeY,0,histsizeY);
  TH3D * h3 = new TH3D("h3", "h3",20,40,60,20,40,60,histsizeZ,0,histsizeZ);


};
