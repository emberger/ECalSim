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

class TROOTAnalysis{

public:
  TROOTAnalysis(TChain* ch);
  ~TROOTAnalysis();

  void plotEvent(Int_t pev);
  void GetCOG(Int_t maxlayer, Int_t cev);


private:

  TTree* EcalTree;
  std::vector<std::pair<Double_t,Double_t>> coglist;



};
