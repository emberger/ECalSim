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

#include "Minimizer.hh"

using namespace ROOT::Minuit2;

class TROOTAnalysis{

public:
  TROOTAnalysis(TChain* ch);
  ~TROOTAnalysis();

  void plotEvent(Int_t pev);
  void CalcCOG(Int_t minlayer, Int_t maxlayer);
  void FitCOGs();
  void PrintFitParams();
  void PrintFitHists();
  //std::vector GetCOGs();


private:
  Int_t nofEntries;
  TTree* EcalTree;
  std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t>> coglist;
  std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t>>> COGCollection;
  std::vector<std::tuple<Double_t, Double_t, Double_t, Double_t>> FitParams;
  std::vector<std::tuple<Double_t, Double_t, Double_t>> ClusteredHits;


};
