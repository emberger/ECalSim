#include "ROOTAnalysis.hh"

TROOTAnalysis::TROOTAnalysis(TChain* ch){
  this->EcalTree=ch->GetTree();
  std::cout << "assigned Tree" << '\n';
  //EcalTree->Print();
}


TROOTAnalysis::~TROOTAnalysis(){}


void TROOTAnalysis::plotEvent(Int_t pev){

  B4ROOTEvent * event = new B4ROOTEvent();
  //std::cout << "Event" << '\n';
  EcalTree->SetBranchAddress("EventBranch", &event);
  //std::cout << "Branch" << '\n';

  Int_t nent=EcalTree->GetEntries();

  if(nent<pev+1){
    std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nent<<"."<<std::endl;
    return;
  }

  TH3D * h = new TH3D("X","ECalEvent",100,50,150,100,50,150,50,0,50);

  //for(Int_t i=0;i<nent; i++){

    EcalTree->GetEntry(pev);
    Int_t pnh=event->NHits();

    for(Int_t j=0;j<pnh;j++){
      h->Fill(event->Hit(j)->X(), event->Hit(j)->Y(),event->Hit(j)->Z(), event->Hit(j)->EnergyDeposit());
    }
//  }
  // const Int_t Number = 3;
  // Double_t Red[Number]    = { 1.00, 0.00, 0.00};
  // Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  // Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
  // Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  // Int_t nb=50;
  // TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  // h->SetContour(nb);
  // h->SetLineWidth(1);
  // h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle("X");
  h->GetYaxis()->SetTitle("Y");
  h->GetZaxis()->SetTitle("Z");

  h->Draw("BOX");
}

void TROOTAnalysis::GetCOG(Int_t maxlayer, Int_t cev){
  B4ROOTEvent * Cevent = new B4ROOTEvent();
  EcalTree->SetBranchAddress("EventBranch", &Cevent);

  EcalTree->GetEntry(cev);

  Int_t cnh=Cevent->NHits();
  Double_t maxDep=0;
  std::pair<Int_t,Int_t> maxTile;

  for(Int_t i=0;i<maxlayer;i++){

    for(Int_t j=0;j<cnh;j++){
      if(Cevent->Hit(j)->Z()==i){
        if(Cevent->Hit(j)->EnergyDeposit()>maxDep){
          maxDep=Cevent->Hit(j)->EnergyDeposit();
          maxTile=std::make_pair(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y());

        }
      }
    }
    std::cout<<maxDep<<" : "<<maxTile.first<<maxTile.second<<std::endl;
    maxDep=0;
    maxTile=std::make_pair(0,0);
  }





}
