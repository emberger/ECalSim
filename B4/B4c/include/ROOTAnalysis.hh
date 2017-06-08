#include "B4ROOTEvent.hh"
#include "TTree.h"
#include "TChain.h"

class ROOTAnalysis{

public:
  ROOTAnalysis(TChain* ch);
  ~ROOTAnalysis();

  void Analyze();


private:

  TTree* EcalTree;



};
