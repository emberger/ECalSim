#include "ROOTAnalysis.hh"

TROOTAnalysis::TROOTAnalysis(TChain* ch){
  this->EcalTree=ch->GetTree();
  std::cout << "assigned Tree" << '\n';
  //EcalTree->Print();
}


TROOTAnalysis::~TROOTAnalysis(){}


void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event

  B4ROOTEvent * event = new B4ROOTEvent();
  //std::cout << "Event" << '\n';
  EcalTree->SetBranchAddress("EventBranch", &event);
  //std::cout << "Branch" << '\n';

  Int_t nent=EcalTree->GetEntries();

  if(nent<pev+1){
    std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nent<<"."<<std::endl;
    return;
  }

  TH3D * h = new TH3D("ECalEvent","ECalEvent",200,0,200,200,0,200,50,0,50);

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

void TROOTAnalysis::GetCOG(Int_t maxlayer, Int_t cev){            //Get vector of (X,Y) tuples containing layerwise center of gravity
  B4ROOTEvent * Cevent = new B4ROOTEvent();
  EcalTree->SetBranchAddress("EventBranch", &Cevent);
  //Int_t nent = EcalTree->GetEntries();

  // EcalTree->GetEntry(cev);
  //
  // Int_t cnh=Cevent->NHits();
  // Double_t maxDep=0;
  //
  //
  //
  // std::tuple<Int_t,Int_t,Double_t, Int_t> seedTiles;
  //
  // for(Int_t i=0;i<maxlayer;i++){
  //
  //   for(Int_t j=0;j<cnh;j++){
  //     if(Cevent->Hit(j)->Z()==i){
  //       std::cout<<"X: "<<Cevent->Hit(j)->X()<<" Y: "<<Cevent->Hit(j)->Y()<<" Z: "<<Cevent->Hit(j)->Z()<<" HitNr: "<<j<<std::endl;
  //       if(Cevent->Hit(j)->EnergyDeposit()>maxDep){
  //         maxDep=Cevent->Hit(j)->EnergyDeposit();
  //         seedTiles=std::make_tuple(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),maxDep, j);
  //
  //       }
  //     }
  //   }
  //   std::cout<<"X "<<std::get<0>(seedTiles)<<" Y: "<<std::get<1>(seedTiles)<<" E: "<<std::get<2>(seedTiles)<<" J: "<<std::get<3>(seedTiles)<<std::endl;
  //   maxDep=0;
  //   seedTiles=std::make_tuple(0,0,0,0);

  TH2D * h1 = new TH2D("h1", "h1", 200,0,200,200,0,200);
  TH2D * h2 = new TH2D("h2", "h2", 200,0,200,200,0,200);

  EcalTree->GetEntry(cev);
  Int_t cnh = Cevent->NHits();
  Double_t integSum=0;
  Double_t integral;

  for(Int_t i=0;i<maxlayer;i++){
    for(Int_t j =0;j<cnh;j++){
      if(Cevent->Hit(j)->Z()==i)
        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());

    }
    integral=h1->Integral();
    Int_t maxbin=h1->GetMaximumBin();
    Double_t maxdep=h1->GetBinContent(maxbin);
    Int_t binx, biny, binz;
    h1->GetBinXYZ(maxbin, binx, biny, binz);


    std::cout<<"Maximum at "<<maxdep<<"MeV"<<" : "<<binx<<":"<<biny<<std::endl;
    std::cout<<"Integral: "<<integral<<std::endl;


    //Do spiral for clustering only if energy in Layer
    if(integral>0){
      Int_t xdir=1;       // intergers for direction
      Int_t ydir=0;
      Int_t buf;
      Int_t stepstodo=1;      // how many steps to go before rotation
      Int_t stepsctr=0;   // number of steps since last rotation

      Int_t currX=binx;  // starting the spiral on the histogram maximum
      Int_t currY=biny;

      Double_t esum=0;
      Double_t curre=0;
      esum+=maxdep;
      h2->SetBinContent(currX, currY, maxdep);
      h1->SetBinContent(binx, biny, 0);
      Int_t ctr=0;

      while(esum<integral*0.8 && currX<binx+20){          //do spiral until desired energyfraction is reached
          currX+=xdir;
          currY+=ydir;
          std::cout<<"X: "<<currX<<"Y: "<<currY<<std::endl;
          curre=h1->GetBinContent(currX, currY);
          std::cout<<"curre"<<curre<<std::endl;
          h2->SetBinContent(currX, currY, curre);
          h1->SetBinContent(currX, currY, 0);
          esum+=curre;
          std::cout<<(esum/integral)*100<<" % of layer energy"<<std::endl;
          stepsctr++;
          if(stepsctr==stepstodo){
            stepsctr=0;
            buf=xdir;         //rotate 90 deg
            xdir= -ydir;
            ydir=buf;

            if(ydir==0){      //incremant steps at every iteration along x
              stepstodo++;
            }
          }
          ctr++;
        }
      }





    integSum += integral;
    h1->Reset();
    h2->Reset();


    std::cout<<integSum<<std::endl;


  }
}
