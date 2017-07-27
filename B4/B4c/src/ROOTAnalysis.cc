#include "ROOTAnalysis.hh"


TROOTAnalysis::TROOTAnalysis(TChain* ch){
  this->EcalTree=ch->GetTree();
  std::cout << "assigned Tree" << '\n';
  nofEntries=EcalTree->GetEntries();
  std::cout<<nofEntries<<std::endl;


    Cevent = new B4ROOTEvent();
    EcalTree->SetBranchAddress("EventBranch", &Cevent);

    // for(Int_t i=0;i<nofEntries;i++){
    //   EcalTree->GetEntry(i);
    //   Int_t hitnr = Cevent->NHits();
    //   for(Int_t j=0;j<hitnr;j++){
    //     Cevent->Hit(j)->SetCoordinates(-495 + Cevent->Hit(j)->X() * 10 , -495 + Cevent->Hit(j)->Y() * 10 , -288.2 + Cevent->Hit(j)->Z() * 11.8 );  //convert Copynumber coordinates
    //   }                                                                                                                                         //to geant4 coordinates
    // }


  //EcalTree->Print();
}


TROOTAnalysis::~TROOTAnalysis(){}


void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event
  TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");

  if(nofEntries<pev+1){
    std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
    return;
  }



  //for(Int_t i=0;i<nent; i++){

    EcalTree->GetEntry(pev);
    Int_t pnh=Cevent->NHits();

    for(Int_t j=0;j<pnh;j++){
      h->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
    }

  h->GetXaxis()->SetTitle("X");
  h->GetYaxis()->SetTitle("Y");
  h->GetZaxis()->SetTitle("Z");

  plotcanvas1->cd();
  h->Draw("BOX");
}


void TROOTAnalysis::findShowercenter(Int_t minevent, Int_t maxevent){



  TH1I * ZHits= new TH1I("ZHits", "ZHits",10000,-295,295);

  for(Int_t events=minevent;events<maxevent;events++){
    EcalTree->GetEntry(events);

    Int_t nhits=Cevent->NHits();
    std::cout<<nhits<<std::endl;

    Int_t maxZ=0;
    for(Int_t i = 0 ; i<nhits ; i++){
      Int_t Z=Cevent->Hit(i)->Z();
      ZHits->Fill(Z);
      std::cout<<"fill"<<std::endl;
      // if(tmp<i){maxZ=Z;}
      // tmp=Z;
    }
    std::cout<<"filled"<<std::endl;
    showerCenters.push_back(ZHits->GetMean());

    std::cout<<ZHits->GetMean()<<std::endl;
    ZHits->Reset();
  }

  for(Int_t j=0;j<maxevent-minevent;j++){
    std::cout<<showerCenters[j]-7<<std::endl;
  }
 //ZHits->Draw("");


}


void TROOTAnalysis::CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

  //variables for clustering
  Int_t xdir=1;       // integers for direction
  Int_t ydir=0;
  Int_t buf;
  Int_t stepstodo=1;      // how many steps to go before rotation
  Int_t stepsctr=0;   // number of steps since last rotation

  Int_t currX=0;  // starting the spiral on the histogram maximum
  Int_t currY=0;

  Double_t esum=0;
  Double_t curre=0;

  Int_t maxbin=0;
  Double_t maxdep=0;
  Int_t binx, biny, binz;

  // Variables for fit
  Double_t xerr=0;
  Double_t yerr=0;
  Double_t cgx=0;
  Double_t cgy=0;
  Double_t cgz=0;
  Double_t Eweight=0;

  for(Int_t eventstodo=minevent; eventstodo<maxevent ;eventstodo++){           //loop over all simulated events

    EcalTree->GetEntry(eventstodo);                   //grab evetn from tree
    Eges = Cevent->GapEnergy();
    Int_t cnh = Cevent->NHits();
    Double_t integral;
    Bool_t foundstart=false;
  for(Int_t i=minlayer;i<maxlayer;i++){               //loop over all layers in event
    for(Int_t j =0;j<cnh ;j++){                       //loop over all hits in laver i
      if(Cevent->Hit(j)->Z()==i)
        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());
         //fill histogram with hits of layer i
    }

    integral=h1->Integral();

    maxbin=h1->GetMaximumBin();
    maxdep=h1->GetBinContent(maxbin);                   // get coordinates of bin with maximum energy deposition
    h1->GetBinXYZ(maxbin, binx, biny, binz);

    //Do spiral for clustering only if energy in Layer

    if(integral!=0){                          //check if layer containes energy

      xdir=1;                                  // integers for direction
      ydir=0;
      stepstodo=1;                             // how many steps to go before rotation
      stepsctr=0;                             // number of steps since last rotation

      currX=binx;                              // starting the spiral on the histogram maximum
      currY=biny;

      esum=0;
      curre=0;
      esum+=maxdep;
      h2->SetBinContent(currX, currY, maxdep);
      auto tp=std::make_tuple(currX, currY, maxdep);
      ClusteredHits.push_back(tp);
      h1->SetBinContent(binx, biny, 0);

      while(esum<integral*0.9 || currX<binx+4){
                                                          //do spiral until desired energyfraction is reached
          currX+=xdir;
          currY+=ydir;
          curre=h1->GetBinContent(currX, currY);

          if(curre>0){                                    //save tile only if it containes energy
            tp=std::make_tuple(currX, currY, curre);
            ClusteredHits.push_back(tp);
            h2->SetBinContent(currX, currY, curre);
          }
          h1->SetBinContent(currX, currY, 0);
          esum+=curre;
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
        }

      Double_t e=h2->Integral();        //get energy contained in cluster

      Eweight=e/Eges;                   //calculate weight

      cgx=h2->GetMean(1);
      xerr=h2->GetMeanError(1);
                                        //extract center of gravity and error
      cgy=h2->GetMean(2);
      yerr=h2->GetMeanError(2);

      if(yerr<0.00001){yerr=1/TMath::Sqrt(12);}
      if(xerr<0.00001){xerr=1/TMath::Sqrt(12);}

      cgz=i;

      er1->Fill(i, xerr/nofEntries);
      er2->Fill(i, yerr/nofEntries);

      h3->Fill(cgx, cgy, cgz);

      auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
      coglist.push_back(cg);

    }         // end of clustering

    ClusteredHits.clear();

    h1->Reset();

    h2->Reset();

  }
  COGCollection.push_back(coglist);
  coglist.clear();

}
  plotCOGs();
}

void TROOTAnalysis::CalcCOG(Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

  //variables for clustering
  Int_t xdir=1;       // integers for direction
  Int_t ydir=0;
  Int_t buf;
  Int_t stepstodo=1;      // how many steps to go before rotation
  Int_t stepsctr=0;   // number of steps since last rotation

  Int_t currX=0;  // starting the spiral on the histogram maximum
  Int_t currY=0;

  Double_t esum=0;
  Double_t curre=0;

  Int_t maxbin=0;
  Double_t maxdep=0;
  Int_t binx, biny, binz;

  // Variables for fit
  Double_t xerr=0;
  Double_t yerr=0;
  Double_t cgx=0;
  Double_t cgy=0;
  Double_t cgz=0;
  Double_t Eweight=0;

  for(Int_t eventstodo=minevent; eventstodo<maxevent ;eventstodo++){           //loop over all simulated events
    std::cout<<"CenterEvent: "<<eventstodo<<std::endl;
    std::cout<<eventstodo-minevent<<std::endl;
    EcalTree->GetEntry(eventstodo);             //grab evetn from tree
    Eges = Cevent->GapEnergy();
    Int_t cnh = Cevent->NHits();
    Double_t integSum=0;
    Double_t integral;
    Bool_t foundstart=false;
    Int_t i = showerCenters[eventstodo-minevent]-7;
    std::cout<<i<<std::endl;
    Int_t collected=0;
    if(i<0){
      i=0;
      collected=4;
    }


  while(collected<7 && i<50){               //loop over layers until some layers after showerstart
    if(i==49 && collected==0){i=0;}
    std::cout<<"i: "<<i<<std::endl;
    for(Int_t j =0;j<cnh ;j++){                       //loop over all hits in laver i
      if(Cevent->Hit(j)->Z()==i)
        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());  //fill histogram with hits of layer i
    }

    integral=h1->Integral();

    maxbin=h1->GetMaximumBin();
    maxdep=h1->GetBinContent(maxbin);                   // get coordinates of bin with maximum energy deposition
    h1->GetBinXYZ(maxbin, binx, biny, binz);

    //Do spiral for clustering only if energy in Layer

    if(integral!=0){                          //check if layer containes energy
      collected++;
      std::cout<<collected<<":"<<i<<std::endl;
      xdir=1;                                  // integers for direction
      ydir=0;
      stepstodo=1;                             // how many steps to go before rotation
      stepsctr=0;                             // number of steps since last rotation

      currX=binx;                              // starting the spiral on the histogram maximum
      currY=biny;

      esum=0;
      curre=0;
      esum+=maxdep;
      h2->SetBinContent(currX, currY, maxdep);
      auto tp=std::make_tuple(currX, currY, maxdep);
      ClusteredHits.push_back(tp);
      h1->SetBinContent(binx, biny, 0);

      while(esum<integral*0.9 || currX<binx+4){
                                                          //do spiral until desired energyfraction is reached
          currX+=xdir;
          currY+=ydir;
          curre=h1->GetBinContent(currX, currY);

          if(curre>0){                                    //save tile only if it containes energy
            tp=std::make_tuple(currX, currY, curre);
            ClusteredHits.push_back(tp);
            h2->SetBinContent(currX, currY, curre);
          }
          h1->SetBinContent(currX, currY, 0);
          esum+=curre;
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
        }

      Double_t e=h2->Integral();        //get energy contained in cluster

      Eweight=e/Eges;                   //calculate weight

      cgx=h2->GetMean(1);
      xerr=h2->GetMeanError(1);
                                        //extract center of gravity and error
      cgy=h2->GetMean(2);
      yerr=h2->GetMeanError(2);

      if(yerr<0.00001){yerr=1/TMath::Sqrt(12);}
      if(xerr<0.00001){xerr=1/TMath::Sqrt(12);}

      cgz=i;

      er1->Fill(i, xerr/nofEntries);
      er2->Fill(i, yerr/nofEntries);

      //h3->Fill(cgx, cgy, cgz);

      auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
      coglist.push_back(cg);

    }         // end of clustering

    integSum += integral;

    ClusteredHits.clear();

    h1->Reset();

    h2->Reset();

    i++;
  }
  COGCollection.push_back(coglist);
  coglist.clear();
}

}

void TROOTAnalysis::CalcCOGwithFit(Int_t minevent, Int_t maxevent){

  TCanvas * c4 = new TCanvas("ReClustering","ReClustering");
  findShowercenter(minevent, maxevent);
  CalcCOG(minevent, maxevent);                          //overloaded CalcCOG function for getting COGs after cluster start
  std::cout<<"center done"<<std::endl;
  FitCOGs(minevent, maxevent);
  std::cout<<"center fit done"<<std::endl;
  //PrintFitHists();
  COGCollection.clear();      // clear the cluster and COG vectors for reclustering
  coglist.clear();
  ClusteredHits.clear();


  //variables for clustering
  Int_t xdir=1;       // integers for direction
  Int_t ydir=0;
  Int_t buf;
  Int_t stepstodo=1;      // how many steps to go before rotation
  Int_t stepsctr=0;   // number of steps since last rotation

  Int_t currX=0;  // starting the spiral on the histogram maximum
  Int_t currY=0;

  Double_t esum=0;
  Double_t curre=0;

  Int_t maxbin=0;
  Double_t maxdep=0;
  Int_t binx, biny, binz;

  // Variables for fit
  Double_t xerr=0;
  Double_t yerr=0;
  Double_t cgx=0;
  Double_t cgy=0;
  Double_t cgz=0;
  Double_t Eweight=0;

  for(Int_t eventstodo=minevent;eventstodo<maxevent;eventstodo++){           //loop over all simulated events
    std::cout<<"ReClustering event: "<<eventstodo<<std::endl;
    EcalTree->GetEntry(eventstodo);             //grab evetn from tree
    Eges = Cevent->GapEnergy();
    //std::cout<<Eges<<std::endl;
    Int_t cnh = Cevent->NHits();
    Double_t integSum=0;
    Double_t integral;


  for(Int_t i=0;i<50;i++){
  //  std::cout<<"layer: "<<i<<std::endl;
    Int_t nofFills=0;                                 //loop over all layers in event
    for(Int_t j =0;j<cnh ;j++){                       //loop over all hits in event i
      if(Cevent->Hit(j)->Z()==i){
        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());  //fill histogram with hits of layer i
        nofFills++;
      }
    }
    //std::cout<<nofFills<<std::endl;

    Double_t MoliereRaduis=48; //cm

    Int_t nofLoops=0;

    while(nofLoops<nofFills){                                        // check if maximum is within 1 MoliereRaduis of Fit

      maxbin=h1->GetMaximumBin();
      maxdep=h1->GetBinContent(maxbin);                              // get coordinates of bin with maximum energy deposition
      h1->GetBinXYZ(maxbin, binx, biny, binz);
      Double_t Xpred=(std::get<0>(FitParams[eventstodo - minevent])*i+std::get<1>(FitParams[eventstodo - minevent]));
      Double_t Ypred=(std::get<2>(FitParams[eventstodo - minevent])*i+std::get<3>(FitParams[eventstodo - minevent]));

    //  std::cout<<"Maxdep: "<<binx<<":"<<biny<<"Pred: "<<Xpred<<":"<<Ypred<<std::endl;

      Double_t d = TMath::Sqrt(    ((-495 + binx*10) - Xpred) * ((-495 + binx*10) - Xpred)  +   ((-495 + biny*10) - Ypred) * ((-495 + biny*10)-Ypred)   );

    //  std::cout<<"d: "<<d<<std::endl;

      if(d>MoliereRaduis){

        h1->SetBinContent(binx, biny, 0);

      }
      else{
    //    std::cout<<"accepted: "<<d<<std::endl;
        break;
      }

      nofLoops++;
    }

    integral=h1->Integral();

    //Do spiral for clustering only if energy in Layer

    if(integral!=0){              //check if layer containes energy

      xdir=1;            // integers for direction
      ydir=0;
      stepstodo=1;       // how many steps to go before rotation
      stepsctr=0;        // number of steps since last rotation

      currX=binx;        // starting the spiral on the histogram maximum
      currY=biny;

      esum=0;
      curre=0;
      esum+=maxdep;
      h2->SetBinContent(currX, currY, maxdep);
      auto tp=std::make_tuple(currX, currY, maxdep);
      ClusteredHits.push_back(tp);
      h1->SetBinContent(binx, biny, 0);

      while(esum<integral*0.9 || currX<binx+9){
                                                          //do spiral until desired energyfraction is reached
          currX+=xdir;
          currY+=ydir;

          curre=h1->GetBinContent(currX, currY);

          if(curre>0){                                    //save tile only if it containes energy
            tp=std::make_tuple(currX, currY, curre);
            ClusteredHits.push_back(tp);
            h2->SetBinContent(currX, currY, curre);
            h1->SetBinContent(currX, currY, 0);
          }


          esum+=curre;
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
        }



      Double_t e=h2->Integral();        //get energy contained in cluster

      Eweight=e/Eges;                   //calculate weight

      cgx=h2->GetMean(1);
      xerr=h2->GetMeanError(1);
                                        //extract center of gravity and error
      cgy=h2->GetMean(2);
      yerr=h2->GetMeanError(2);

      if(yerr<0.00001){yerr=1/TMath::Sqrt(12);}
      if(xerr<0.00001){xerr=1/TMath::Sqrt(12);}

      cgz=i;

      er1->Fill(i, xerr/nofEntries);
      er2->Fill(i, yerr/nofEntries);



      //std::cout<<"Center of gravity: "<<cgx<<" "<<cgy<<" "<<cgz<<std::endl;
      h3->Fill(cgx, cgy, cgz);

      auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
      coglist.push_back(cg);

    }         // end of clustering


    integSum += integral;



    ClusteredHits.clear();

    h1->Reset();

    h2->Reset();

    Eweight=0;




  }                   // loop over layers

  COGCollection.push_back(coglist);
  coglist.clear();

}                     // loop over events
  FitParams.clear();
  FitCOGs(minevent, maxevent);
  PrintFitHists(minevent, maxevent);
  plotCOGs();
  //std::cout<<breakctr<<std::endl;
}



void TROOTAnalysis::plotCOGs(){

  h3->GetXaxis()->SetTitle("X");
  h3->GetYaxis()->SetTitle("Y");
  h3->GetZaxis()->SetTitle("Z");
  h3->SetMarkerStyle(4);
  c2->cd();
  h3->Draw("");
  //h3->Reset();

}



void TROOTAnalysis::FitCOGs( Int_t minevent, Int_t maxevent){
  //histograms for correlation
  TCanvas * corr1 = new TCanvas("Correlations");
  corr1->Divide(3,2,0.01,0.01);

  TH1D * co1 = new TH1D("MxTx correllation", "MxTx correllation",100, -1,1 );
  co1->GetXaxis()->SetTitle("corellation of X slope and X intercept");
  TH1D * co2 = new TH1D("MxMy correllation", "MxMy correllation",100, -1,1 );
  co2->GetXaxis()->SetTitle("corellation of X slope and Y slope");
  TH1D * co3 = new TH1D("MxTy correllation", "MxTy correllation",100, -1,1 );
  co3->GetXaxis()->SetTitle("correlation of X slope and Y intercept");
  TH1D * co4 = new TH1D("TxMy correllation", "TxMy correllation",100, -1,1 );
  co4->GetXaxis()->SetTitle("correlation of X intercept and Y slope");
  TH1D * co5 = new TH1D("TxTy correllation", "TxTy correllation",100, -1,1 );
  co5->GetXaxis()->SetTitle("correlation of X intercapt an Y intercept");
  TH1D * co6 = new TH1D("MyTy correllation", "MyTy correllation",100, -1,1 );
  co6->GetXaxis()->SetTitle("correlation of Y slope and Y intercept");

  TCanvas * dist1= new TCanvas("2m distance");
  dist1->Divide(2,1,0.01,0.01);
  TH1D * dx = new TH1D("Front X", "Front X", 2000, -50,150);
  dx->GetXaxis()->SetTitle("X[cm]");
  TH1D * dy = new TH1D("Front Y", "Front Y", 2000, -50,150);
  dy->GetXaxis()->SetTitle("Y[cm]");


    //transform from copynumber coordinates to geant4 coordinates
    std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t>> Transcoglist;
    std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t>>> TransfomedCOGs;

    for(Int_t i = 0;i<nofEntries;i++){
      for(Int_t j = 0;j<COGCollection[i].size();j++){

        auto tc = std::make_tuple((-495 + std::get<0>(COGCollection[i][j]) * 10) ,
                                  (-495 + std::get<1>(COGCollection[i][j]) * 10) ,
                                  -288.2 + std::get<2>(COGCollection[i][j]) * 11.8,
                                  std::get<3>(COGCollection[i][j]),
                                  std::get<4>(COGCollection[i][j]),
                                  std::get<5>(COGCollection[i][j])  );
        Transcoglist.push_back(tc);

      }
      TransfomedCOGs.push_back(Transcoglist);
      Transcoglist.clear();
    }

    for(Int_t i = 0; i<nofEntries;i++){
      for(Int_t j =0;j<TransfomedCOGs[i].size();j++){
        std::cout<<std::get<0>(TransfomedCOGs[i][j])<<":"<<std::get<1>(TransfomedCOGs[i][j])<<":"<<std::get<2>(TransfomedCOGs[i][j])<<"-----"
        <<std::get<0>(COGCollection[i][j])<<":"<<std::get<1>(COGCollection[i][j])<<":"<<std::get<2>(COGCollection[i][j])<<std::endl;;
      }
    }


    Fcn myfcn;
    myfcn.SetCOGs(TransfomedCOGs);
    //myfcn.PrintCOGs();
    MnUserParameters upar;
    double error_minimizer_parameters = 1e-4;

    upar.Add("mx", 0., error_minimizer_parameters);
    upar.Add("tx", 1., error_minimizer_parameters);
    upar.Add("my", 0., error_minimizer_parameters);
    upar.Add("ty", 1., error_minimizer_parameters);

    cout << upar << endl;

    MnMigrad migrad(myfcn, upar, 2);



    for(int events = 0; events < maxevent-minevent; events++){

        upar.SetValue("mx", 0.);
        upar.SetValue("tx", 0.);
        upar.SetValue("my", 0.);
        upar.SetValue("ty", 0.);

        myfcn.SetCurrentEvent(events);

        std::cout<<"event: "<<events<<std::endl;
        FunctionMinimum min = migrad();

        //std::cout<<min<<std::endl;


        MnUserCovariance cov = min.UserCovariance();

        //std::cout<<cov(0 ,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<cov(3,3)<<std::endl;


        MnUserParameterState userParameterState = min.UserState();

        auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                        userParameterState.Value("my"), userParameterState.Value("ty"));
        FitParams.push_back(tp);

        dx->Fill(std::get<0>(tp)*(-295)+std::get<1>(tp));
        dy->Fill(std::get<2>(tp)*(-295)+std::get<3>(tp));

        Double_t covar[4][4];
        Double_t error[4];

        error[0]= userParameterState.Error("mx");
        error[1]= userParameterState.Error("tx");
        error[2]= userParameterState.Error("my");
        error[3]= userParameterState.Error("ty");

        //Fill correlation Histograms
        if(COGCollection[events].size() > 1){


          for(Int_t row=0;row<4;row++){
            for(Int_t col=0;col<4;col++){
              covar[row][col]=cov(row,col);

            }
          }

          co1->Fill(covar[0][1]/(error[0]*error[1]));
          co2->Fill(covar[0][2]/(error[0]*error[2]));
          co3->Fill(covar[0][3]/(error[0]*error[3]));
          co4->Fill(covar[1][2]/(error[1]*error[2]));
          co5->Fill(covar[1][3]/(error[1]*error[3]));
          co6->Fill(covar[2][3]/(error[2]*error[3]));
        }
        // auto tp=std::make_tuple(userParameterState.Value("a1"),userParameterState.Value("a2"),userParameterState.Value("a3"),
        //                         userParameterState.Value("x1"),userParameterState.Value("x2"),userParameterState.Value("x3"));
        //
        // FitParams.push_back(tp);
    }

    dist1->cd(1);
    dx->Draw("");

    dist1->cd(2);
    dy->Draw("");

    corr1->cd(1);
    co1->Draw("");

    corr1->cd(2);
    co2->Draw("");

    corr1->cd(3);
    co3->Draw("");

    corr1->cd(4);
    co4->Draw("");

    corr1->cd(5);
    co5->Draw("");

    corr1->cd(6);
    co6->Draw("");

}


void TROOTAnalysis::PrintFitHists(Int_t minevent, Int_t maxevent){
  Int_t ent=maxevent-minevent;
  TCanvas * c1 = new TCanvas("Fit", "Fit");
  c1->Divide(2,3,0.01, 0.01);



  for(Int_t i=0;i<ent;i++){
      fitX->Fill(std::get<1>(FitParams[i]));
  }

  for(Int_t i=0;i<ent;i++){
      fitY->Fill(std::get<3>(FitParams[i]));
  }

  for(Int_t i=0;i<ent;i++){
      fitSX->Fill(std::get<0>(FitParams[i]));
  }

  for(Int_t i=0;i<ent;i++){
      fitSY->Fill(std::get<2>(FitParams[i]));
  }

  c1->cd(1);
  fitX->GetXaxis()->SetTitle("X");
  fitX->GetYaxis()->SetTitle("#");
  fitX->Draw();

  c1->cd(2);
  fitY->GetXaxis()->SetTitle("Y");
  fitY->GetYaxis()->SetTitle("#");
  fitY->Draw();

  c1->cd(3);
  fitSX->GetXaxis()->SetTitle("X Slope");
  fitSX->GetYaxis()->SetTitle("#");
  fitSX->Draw();

  c1->cd(4);
  fitSY->GetXaxis()->SetTitle("Y Slope");
  fitSY->GetYaxis()->SetTitle("#");
  fitSY->Draw();

  c1->cd(5);
  er1->Draw();

  c1->cd(6);
  er2->Draw();

//  er1->Reset();
//  er2->Reset();
}

void TROOTAnalysis::CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){



  TCanvas * c3 = new TCanvas("CleanCOGs", "CleanCOGs");
  TH3D * clCOG = new TH3D("COGs after cleaning", "COGs after celaning",20,40,60,20,40,60,50,0,50);

  Double_t MoliereRaduis=4.87;   // in cm

  CalcCOG(0,10, minevent-1, maxevent-1);
  FitCOGs(minevent-1, maxevent-1);

  COGCollection.clear();
  coglist.clear();
  ClusteredHits.clear();

  //PrintFitHists();
  CalcCOG(minlayer-1,maxlayer-1,minevent-1, maxevent-1 );

  Int_t erasecounter=0;
  Int_t totalcounter=0;
  for(int i=0; i<nofEntries;i++){                                       //Delete the COGs, outside of one MoliereRadius
    for(int z=0; z<COGCollection[i].size(); z++){

      Double_t distx=(std::get<0>(FitParams[i])*z+std::get<1>(FitParams[i]))-std::get<0>(COGCollection[i][z]);
      Double_t disty=(std::get<2>(FitParams[i])*z+std::get<3>(FitParams[i]))-std::get<1>(COGCollection[i][z]);
      totalcounter++;
      if(TMath::Sqrt(distx*distx)>MoliereRaduis || TMath::Sqrt(disty*disty)>MoliereRaduis){
        COGCollection[i].erase(COGCollection[i].begin()+z);
        //std::cout<<"X distance: "<<distx<<" Y distance: "<<disty<<std::endl;
        erasecounter++;
      }

    }

  }

  for(int i=0; i<nofEntries;i++){                                       //Delete the COGs, outside of one MoliereRadius
    for(int z=0; z<COGCollection[i].size(); z++){
      clCOG->Fill(std::get<0>(COGCollection[i][z]), std::get<1>(COGCollection[i][z]), std::get<2>(COGCollection[i][z]));
    }
  }

  clCOG->GetXaxis()->SetTitle("X");
  clCOG->GetYaxis()->SetTitle("Y");
  clCOG->GetZaxis()->SetTitle("Z");
  clCOG->SetMarkerStyle(4);
  c3->cd();
  clCOG->Draw("");

  FitParams.clear();


  FitCOGs(minevent, maxevent);                //refit the leftover COGs
  PrintFitHists(minevent, maxevent);
  plotCOGs();

  std::cout<<"Discarded COGs: "<<erasecounter<<" Total COGs:"<<totalcounter<<std::endl;




}
