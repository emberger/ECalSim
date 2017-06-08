#include "ROOTAnalysis.hh"



ROOTAnalysis::ROOTAnalysis(TChain* ch){
  EcalTree=ch->GetTree();

}

ROOTAnalysis::~ROOTAnalysis(){}


void ROOTAnalysis::Analyze(){

  B4ROOTEvent * event = new B4ROOTEvent();

  EcalTree->SetBranchAddress("Event", &event);



}
