#include "ROOTAnalysis.hh"


TROOTAnalysis::TROOTAnalysis(TChain* ch){
  this->EcalTree=ch->GetTree();
  std::cout << "assigned Tree" << '\n';
  nofEntries=EcalTree->GetEntries();
  //EcalTree->Print();
}


TROOTAnalysis::~TROOTAnalysis(){}


void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event

  B4ROOTEvent * event = new B4ROOTEvent();
  //std::cout << "Event" << '\n';
  EcalTree->SetBranchAddress("EventBranch", &event);
  //std::cout << "Branch" << '\n';

  if(nofEntries<pev+1){
    std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
    return;
  }

  TH3D * h = new TH3D("ECalEvent","ECalEvent",500,0,500,500,0,500,50,0,50);

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

void TROOTAnalysis::CalcCOG(Int_t minlayer, Int_t maxlayer){            //Get vector of (X,Y) tuples containing layerwise center of gravity
  B4ROOTEvent * Cevent = new B4ROOTEvent();
  EcalTree->SetBranchAddress("EventBranch", &Cevent);


  Double_t histsizeX=100;
  Double_t histsizeY=100;
  Double_t histsizeZ=50;


  TH2D * h1 = new TH2D("h1", "h1", histsizeX,0,histsizeX,histsizeY,0,histsizeY);
  TH2D * h2 = new TH2D("h2", "h2", histsizeX,0,histsizeX,histsizeY,0,histsizeY);
  TH3D * h3 = new TH3D("h3", "h3",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
  //TGraph2D * g = new TGraph2D();


  for(Int_t eventstodo=0;eventstodo<nofEntries;eventstodo++){
    EcalTree->GetEntry(eventstodo);
    Int_t cnh = Cevent->NHits();
    Double_t integSum=0;
    Double_t integral;
    Int_t graphctr=0;

  for(Int_t i=minlayer;i<maxlayer;i++){
    for(Int_t j =0;j<cnh ;j++){
      if(Cevent->Hit(j)->Z()==i)
        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());

    }
    integral=h1->Integral();

    Int_t maxbin=h1->GetMaximumBin();
    Double_t maxdep=h1->GetBinContent(maxbin);
    Int_t binx, biny, binz;
    h1->GetBinXYZ(maxbin, binx, biny, binz);

    // std::cout<<"Maximum at "<<maxdep<<"MeV"<<" : "<<binx<<":"<<biny<<std::endl;
    //std::cout<<"Integral: "<<integral<<std::endl;


    //Do spiral for clustering only if energy in Layer
    if(integral!=0 && maxdep > 0){


      Int_t xdir=1;       // integers for direction
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
      auto tp=std::make_tuple(currX, currY, maxdep);
      ClusteredHits.push_back(tp);
      h1->SetBinContent(binx, biny, 0);
      Int_t ctr=0;

      while(esum<integral*0.8 && currX<binx+12){          //do spiral until desired energyfraction is reached
          currX+=xdir;
          currY+=ydir;
          //std::cout<<"X: "<<currX<<"Y: "<<currY<<std::endl;
          curre=h1->GetBinContent(currX, currY);
          //std::cout<<"curre"<<curre<<std::endl;
          if(curre>0){
            tp=std::make_tuple(currX, currY, curre);
            ClusteredHits.push_back(tp);
            h2->SetBinContent(currX, currY, curre);
          }
          h1->SetBinContent(currX, currY, 0);
          esum+=curre;
          //std::cout<<(esum/integral)*100<<" % of layer energy"<<std::endl;
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

      // for(Int_t foo=0;foo<ClusteredHits.size();foo++){
      //   std::cout<<"X: "<<std::get<0>(ClusteredHits[foo])<<" Y: "<<std::get<1>(ClusteredHits[foo])<<" E: "<<std::get<2>(ClusteredHits[foo])<<std::endl;
      // }



      Double_t integral2= h2->Integral();

      Double_t cgx=h2->GetMean(1);
      //Double_t xerr=h2->GetMeanError(1);

      Double_t cgy=h2->GetMean(2);
      //Double_t yerr=h2->GetMeanError(2);

      //Double_t cgx=0;
      //Double_t cgy=0;
      Double_t cgz=i;

      Double_t xerr=0;
      Double_t yerr=0;
      Double_t counter=0;
      for(Int_t er=0; er<ClusteredHits.size();er++){
        //cgx += (std::get<0>(ClusteredHits[er])*std::get<2>(ClusteredHits[er]));
        //cgy += (std::get<1>(ClusteredHits[er])*std::get<2>(ClusteredHits[er]));

        xerr+= std::get<2>(ClusteredHits[er])*std::get<2>(ClusteredHits[er]);
        yerr+= std::get<2>(ClusteredHits[er])*std::get<2>(ClusteredHits[er]);
        counter++;
      }



      Double_t tmp1=TMath::Sqrt(xerr)/(TMath::Sqrt(12.0)*integral2);
      Double_t tmp2=TMath::Sqrt(yerr)/(TMath::Sqrt(12.0)*integral2);

      xerr=tmp1;
      yerr=tmp2;


      // std::cout<<"-------"<<std::endl;
      // std::cout<<integral2<<" hits: "<<ClusteredHits.size()<<"Counter: "<<counter<<"energy in first entry"<<std::get<2>(ClusteredHits[0])<<std::endl;
      // std::cout<<xerr<<" : "<<yerr<<std::endl;
      //
      // std::cout<<"Layer "<<i<<"done"<<std::endl;
      // std::cout<<"-------"<<std::endl;
      //cgx=cgx/integral2;
      //cgy=cgy/integral2;


      //std::cout<<"Center of gravity: "<<cgx<<" "<<cgy<<" "<<cgz<<std::endl;
      h3->Fill(cgx, cgy, cgz);
      //g->SetPoint(graphctr, cgx, cgy, cgz);
      if(integral2 != 0){
        auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr);
        coglist.push_back(cg);
      }
    }


    integSum += integral;
    graphctr++;
    ClusteredHits.clear();
    std::vector<std::tuple<Double_t, Double_t, Double_t>>().swap(ClusteredHits);
    h1->Reset();
    h2->Reset();





  }
  //std::cout<<"Event done"<<std::endl;
  //std::cout<<integSum<<std::endl;
  COGCollection.push_back(coglist);
  coglist.clear();
   h3->GetXaxis()->SetTitle("X");
   h3->GetYaxis()->SetTitle("Y");
   h3->GetZaxis()->SetTitle("Z");
   h3->SetMarkerStyle(4);
   //h3->Draw("");
   //h3->Reset();
  //g->Draw("");

}
// for(Int_t p=0;p<nofEntries;p++){
//   for(Int_t q=0;q<COGCollection[p].size();q++){
// std::cout<<std::get<0>(COGCollection[p][q])<<":"<<std::get<1>(COGCollection[p][q])<<":"<<std::get<2>(COGCollection[p][q])<<std::endl;
// std::cout<<"Errorx: "<<std::get<3>(COGCollection[p][q])<<"Errory: "<<std::get<4>(COGCollection[p][q])<<std::endl;
//
// }
// std::cout<<"--------------------"<<std::endl;
// }
}

void TROOTAnalysis::FitCOGs(){


    Fcn myfcn;
    myfcn.SetCOGs(COGCollection);
    //myfcn.PrintCOGs();
    MnUserParameters upar;
    double error_minimizer_parameters = 1e-4;

    upar.Add("mx", 0., error_minimizer_parameters);
    upar.Add("tx", 1., error_minimizer_parameters);
    upar.Add("my", 0., error_minimizer_parameters);
    upar.Add("ty", 1., error_minimizer_parameters);

    cout << upar << endl;

    MnMigrad migrad(myfcn, upar, 2);

    for(int events = 0; events <nofEntries; events++){

        upar.SetValue("mx", 0.);
        upar.SetValue("tx", 0.);
        upar.SetValue("my", 0.);
        upar.SetValue("ty", 0.);

        myfcn.SetCurrentEvent(events);


        FunctionMinimum min = migrad();



        MnUserParameterState userParameterState = min.UserState();

        auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                        userParameterState.Value("my"), userParameterState.Value("ty"));
        FitParams.push_back(tp);



    }
}

void TROOTAnalysis::PrintFitParams(){

  for(Int_t p=0;p<nofEntries;p++){
    std::cout<<"Event 1 ------------------------"<<std::endl;
    std::cout<<"mx: "<<std::get<0>(FitParams[p])<<" tx: "<<std::get<1>(FitParams[p])<<" my: "<<std::get<2>(FitParams[p])<<" ty: "<<std::get<3>(FitParams[p])<<std::endl;

  }



}

void TROOTAnalysis::PrintFitHists(){

  TH1D * fit1 = new TH1D("X Fit","X Fit", 1000,55,65);
  for(Int_t i=0;i<nofEntries;i++){
      fit1->Fill(std::get<1>(FitParams[i]));

  }


  fit1->GetXaxis()->SetTitle("X");
  fit1->GetYaxis()->SetTitle("#");
  fit1->Draw();


}





// std::vector TROOTAnalysis::GetCOGs(){return coglist;}
