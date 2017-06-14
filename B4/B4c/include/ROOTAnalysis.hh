#include "B4ROOTEvent.hh"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TColor.h"
#include "TH3D.h"

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
