#include "TROOTAnalysis.hh"




TROOTAnalysis::TROOTAnalysis(TChain* ch){
        this->EcalTree=ch->GetTree();
        std::cout << "assigned Tree" << '\n';
        nofEntries=EcalTree->GetEntries();
        std::cout<<nofEntries<<std::endl;
        FitParamsPions.resize(nofEntries);

        Cevent = new B4ROOTEvent();
        EcalTree->SetBranchAddress("EventBranch", &Cevent);

        // for(Int_t i=0; i<nofEntries; i++) {                                              //energySmearing
        //         EcalTree->GetEntry(i);
        //         Int_t hitnr = Cevent->NHits();
        //         for(Int_t j=0; j<hitnr; j++) {
        //                 std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
        //         }
        // }


        //EcalTree->Print();
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

TROOTAnalysis::~TROOTAnalysis(){
}


void TROOTAnalysis::plotEvent(Int_t pev){     //plot 3DHisto of selected event
        TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");

        if(nofEntries<pev+1) {
                std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
                return;
        }



        //for(Int_t i=0;i<nent; i++){

        EcalTree->GetEntry(pev);
        Int_t pnh=Cevent->NHits();

        for(Int_t j=0; j<pnh; j++) {
                h->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
        }

        h->GetXaxis()->SetTitle("X");
        h->GetYaxis()->SetTitle("Y");
        h->GetZaxis()->SetTitle("Z");

        plotcanvas1->cd();
        h->Draw("BOX");

        std::string Path="/home/iwsatlas1/emberger/FuerFrank/PhotonEventDisplay/Event";
        std::string nr = std::to_string(pev);

        std::string extension = ".C";
        std::string PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        extension = ".pdf";
        PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());
        h->Reset();
}

void TROOTAnalysis::plotEventPion(Int_t pev){     //plot 3DHisto of selected event
        TCanvas * plotcanvas1 = new TCanvas("eventplotter", "eventplotter");
        if(nofEntries<pev+1) {
                std::cout<<"Tree has less than "<<pev+1<<" events, it has "<<nofEntries<<"."<<std::endl;
                return;
        }



        //for(Int_t i=0;i<nent; i++){

        EcalTree->GetEntry(pev);
        Int_t pnh=Cevent->NHits();

        for(Int_t j=0; j<pnh; j++) {
                if(Cevent->Hit(j)->PhotNr()==1) {
                        hA->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
                }
        }

        //TColor * col1= new TColor(1,0,0,"red");

        hA->GetXaxis()->SetTitle("X");
        hA->GetYaxis()->SetTitle("Y");
        hA->GetZaxis()->SetTitle("Z");
        hA->SetMarkerColor(kRed);
        hA->SetMarkerStyle(20);
        hA->SetMarkerSize(0.3);
        hA->Draw("");

        for(Int_t j=0; j<pnh; j++) {
                if(Cevent->Hit(j)->PhotNr()==2) {
                        hB->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(),Cevent->Hit(j)->Z(), Cevent->Hit(j)->EnergyDeposit());
                }
        }

        //TColor * col2 = new TColor(0,0,1,"blue");
        hB->SetMarkerColor(kBlue);
        hB->SetMarkerStyle(20);
        hB->SetMarkerSize(0.3);
        hB->GetXaxis()->SetTitle("X");
        hB->GetYaxis()->SetTitle("Y");
        hB->GetZaxis()->SetTitle("Z");
        hB->Draw("SAME");

        std::string Path="/home/iwsatlas1/emberger/FuerFrank/PionEventDisplay/Event";
        std::string nr = std::to_string(pev);

        std::string extension = ".C";
        std::string PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        // extension = ".png";
        // PlotPath=Path+nr+extension;
        // plotcanvas1->Print(PlotPath.c_str());

        extension = ".pdf";
        PlotPath=Path+nr+extension;
        plotcanvas1->Print(PlotPath.c_str());

        hB->Reset();
        hA->Reset();
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotRMSx(){
        TCanvas * rms=new TCanvas("rmsx_over_E");


        // Double_t y[10]={0.5298,0.3507,0.3624,0.2524,0.222,0.1861,0.228,0.1751,0.1425,0.1481}; // 0.0 slope
        // Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        //
        Double_t y[10]={0.5398,0.3428,0.3087,0.2319,0.2117,0.1933,0.2243,0.1716,0.1411,0.1588};                                                                      //0.4slope
        Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        TGraph *gRMS1=new TGraph(10, x, y );
        gRMS1->Fit(rms_x);
        gStyle->SetOptFit(11111111);
        gRMS1->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS1->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS1->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(0);
        gRMS1->Draw("A*");
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PrintERes(){

        TCanvas * res = new TCanvas("Energy Resolution", "ERes",2000,1000);
        TCanvas * gap = new TCanvas("GapEnergy", "GapEnergy");

        Double_t y[12]={4.55/19.25, 5.942/38.47, 9.056/77.51, 12.06/115.9, 11.59/154.2, 13.71/192.9, 17.98/268.8, 20.25/384, 53.25/956.2, 69.23/1905, 94.6/2836, 146.7/3779};
        Double_t x[12]={0.05,0.100,0.200,0.300,0.400,0.500,0.700,1.000,2.500,5.000,7.500,10.000};

        TGraph * re1 = new TGraph(12, x, y);
        re1->SetTitle("Energyresolution");
        TF1 * EnergyRes = new TF1("EnergyRes", "sqrt(([0] * [0] / x) + ([1]*[1])  + ([2]*[2]/(x*x)))");

        re1->Fit(EnergyRes);



        re1->GetXaxis()->SetTitle("Energy[GeV]");
        re1->GetYaxis()->SetTitle("#frac{#sigma}{E}");
        re1->GetYaxis()->SetTitleOffset(1.3);
        //re1->GetYaxis()->LabelsOption("v");
        //re1->SetMarkerStyle(23);
        gStyle->SetOptFit();
        res->cd(0);
        re1->Draw("A*");
        TImage * img1 = TImage::Create();
        img1->FromPad(res);
        //img1->Scale(1000, 1000);
        img1->WriteImage("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/GammaEnergyandSlopeScan_Analysis/ERes.png");
        delete img1;

        Double_t Egun[12]={50,100,200,300,400,500,700,1000,2500,5000,7500,10000};
        Double_t Egap[12]={19.25,38.47,77.51,115.9,154.2,192.9,268.8,384,956.2,1905,2836,3779};

        Double_t errx[12]={0,0,0,0,0,0,0,0,0,0,0,0};
        Double_t erry[12]={4.55, 5.942, 9.056, 12.06, 11.59, 13.71, 17.98, 20.25, 53.25, 69.23, 94.6, 146.7};

        TF1 * Efficiency = new TF1("Efficiency", "[0]*x+[1]");

        TGraphErrors * gap1 = new TGraphErrors(12, Egun, Egap, errx, erry );
        gap1->SetTitle("Efficiency");
        gap1->Fit(Efficiency);
        gap1->GetXaxis()->SetTitle("Gun Energy[MeV]");
        gap1->GetYaxis()->SetTitle("Gap Energy[MeV]");
        gap1->GetYaxis()->SetTitleOffset(1.3);
        gStyle->SetOptFit();

        gap->cd(0);
        gap1->Draw("A*");
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::PlotProjection(Int_t distance){

        TCanvas * dist1= new TCanvas("Projection", "Projection", 1500,1000);
        TCanvas * dist2 = new TCanvas("ProjectionXYC", "ProjectionXYC", 1000,1000);


        dist1->Divide(2,1,0.01,0.01);

        TH1D * dx = new TH1D("Projection X", "Projection X", 1000, -600,600);
        dx->GetXaxis()->SetTitle("X[mm]");
        TH1D * dy = new TH1D("Projection Y", "Projection Y", 1000, -600,600);
        dy->GetXaxis()->SetTitle("Y[mm]");

        TH2D * cont1 = new TH2D("ProjectionXY", "ProjectionXY", 2000, -1000,1000,2000,-1000,1000);
        cont1->GetXaxis()->SetTitle("X[mm]");
        cont1->GetYaxis()->SetTitle("Y[mm]");
        std::cout<<"Created Hists"<<std::endl;

        for(Int_t i=0; i<nofEntries; i++) {
                Double_t proX=std::get<0>(FitParamsGamma[i])*(-295-distance)+std::get<1>(FitParamsGamma[i]);
                Double_t proY=std::get<2>(FitParamsGamma[i])*(-295-distance)+std::get<3>(FitParamsGamma[i]);

                dx->Fill(proX);
                dy->Fill(proY);
                cont1->Fill(proX,proY);
                std::cout<<"filled: "<<i<<std::endl;
        }

        gStyle->SetStatY(0.9);
        // Set y-position (fraction of pad size)
        gStyle->SetStatX(0.9);
        // Set x-position (fraction of pad size)

        dist1->cd(1);
        dx->Draw("");

        dist1->cd(2);
        dy->Draw("");

        dist2->cd(0);
        gStyle->SetStatY(0.9);
        // Set y-position (fraction of pad size)
        gStyle->SetStatX(0.9);
        // Set x-position (fraction of pad size)
        cont1->Draw("colz");


        if(pathset) {

                // std::string imgname1= "Projection100mm.png";
                // std::string imgpath1 = savepath + '/' + imgname1;

                std::string imgname2= "ProjectionColour100mm";
                std::string imgpath2c = savepath + '/' + imgname2;


                std::string imgpath2=imgpath2c + ".C";
                dist2->Print(imgpath2.c_str());


                imgpath2=imgpath2c + ".pdf";
                dist2->Print(imgpath2.c_str());




                // TImage * img3 = TImage::Create();
                // TImage * img4 = TImage::Create();
                //
                // img3->FromPad(dist1);
                // img4->FromPad(dist2);
                //
                // img3->WriteImage(imgpath1.c_str());
                // img4->WriteImage(imgpath2.c_str());
                // delete img3;
                // delete img4;
        }

}


//--------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::plotCOGs(){

        h3->GetXaxis()->SetTitle("X");
        h3->GetYaxis()->SetTitle("Y");
        h3->GetZaxis()->SetTitle("Z");
        h3->SetMarkerStyle(4);
        c2->cd();
        h3->Draw("");
        //h3->Reset();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------




void TROOTAnalysis::findShowercenter(Int_t minevent, Int_t maxevent){



        TH1I * ZHits= new TH1I("ZHits", "ZHits",10000,-295,295);

        for(Int_t events=minevent; events<maxevent; events++) {
                EcalTree->GetEntry(events);

                Int_t nhits=Cevent->NHits();
                std::cout<<nhits<<std::endl;

                Int_t maxZ=0;

                for(Int_t i = 0; i<nhits; i++) {
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

        for(Int_t j=0; j<maxevent-minevent; j++) {
                std::cout<<showerCenters[j]-7<<std::endl;
        }
        //ZHits->Draw("");


}

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::CalcCOGPion(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent,  Int_t photonNR){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

        //variables for clustering
        Int_t xdir=1; // integers for direction
        Int_t ydir=0;
        Int_t buf;
        Int_t stepstodo=1; // how many steps to go before rotation
        Int_t stepsctr=0; // number of steps since last rotation

        Int_t currX=0; // starting the spiral on the histogram maximum
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

        for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) { //loop over all simulated events

                EcalTree->GetEntry(eventstodo); //grab evetn from tree
                Eges = Cevent->GapEnergy();
                Int_t cnh = Cevent->NHits();
                Double_t integral;
                for(Int_t i=minlayer; i<maxlayer; i++) { //loop over all layers in event
                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i
                                //std::cout<<Cevent->Hit(j)->PhotNr()<<std::endl;
                                if(Cevent->Hit(j)->Z()==i && Cevent->Hit(j)->PhotNr()==photonNR) {
                                        //std::cout<<"fill: "<<Cevent->Hit(j)->PhotNr()<<"from photon mode:"<<photonNR<<std::endl;
                                        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());
                                }//fill histogram with hits of layer i
                        }
                        //std::cout<<photonNR<<" : "<<std::endl;
                        integral=h1->Integral();

                        maxbin=h1->GetMaximumBin();
                        maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
                        h1->GetBinXYZ(maxbin, binx, biny, binz);

                        //Do spiral for clustering only if energy in Layer

                        if(integral!=0) { //check if layer containes energy

                                xdir=1;  // integers for direction
                                ydir=0;
                                stepstodo=1; // how many steps to go before rotation
                                stepsctr=0; // number of steps since last rotation

                                currX=binx; // starting the spiral on the histogram maximum
                                currY=biny;

                                esum=0;
                                curre=0;
                                esum+=maxdep;
                                h2->SetBinContent(currX, currY, maxdep);
                                auto tp=std::make_tuple(currX, currY, maxdep);
                                ClusteredHits.push_back(tp);
                                h1->SetBinContent(binx, biny, 0);

                                while(esum<integral*0.9 || currX<binx+4) {
                                        //do spiral until desired energyfraction is reached
                                        currX+=xdir;
                                        currY+=ydir;
                                        curre=h1->GetBinContent(currX, currY);

                                        if(curre>0) { //save tile only if it containes energy
                                                tp=std::make_tuple(currX, currY, curre);
                                                ClusteredHits.push_back(tp);
                                                h2->SetBinContent(currX, currY, curre);
                                        }
                                        h1->SetBinContent(currX, currY, 0);
                                        esum+=curre;
                                        stepsctr++;
                                        if(stepsctr==stepstodo) {
                                                stepsctr=0;
                                                buf=xdir; //rotate 90 deg
                                                xdir= -ydir;
                                                ydir=buf;

                                                if(ydir==0) { //incremant steps at every iteration along x
                                                        stepstodo++;
                                                }
                                        }
                                }

                                Double_t e=h2->Integral(); //get energy contained in cluster

                                Eweight=e/Eges; //calculate weight

                                cgx=h2->GetMean(1);
                                xerr=h2->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h2->GetMean(2);
                                yerr=h2->GetMeanError(2);

                                if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                                if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                                cgz=i;

                                er1->Fill(i, xerr/nofEntries);
                                er2->Fill(i, yerr/nofEntries);

                                h3->Fill(cgx, cgy, cgz);

                                auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
                                coglist.push_back(cg);

                        } // end of clustering

                        ClusteredHits.clear();

                        h1->Reset();

                        h2->Reset();

                }
                COGCollection.push_back(coglist);
                coglist.clear();

        }
        plotCOGs();
}


//-------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::AnalyzePions(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){

        Bool_t pion=true;
        CalcCOGPion(minlayer, maxlayer, minevent, maxevent, 1);
        plotCOGs();
        FitCOGsPion(minevent, maxevent, 10, pion, 1);

        COGCollection.clear();
        coglist.clear();

        CalcCOGPion(minlayer, maxlayer, minevent, maxevent, 2);
        FitCOGsPion(minevent, maxevent, 10, pion, 2);
        //
        // for(Int_t i=0; i<nofEntries; i++) {
        //         std::cout<<"Photon1: "<<std::endl;
        //         std::cout<<std::get<0>(FitParamsPions[i][0])<<"*z + "<<std::get<1>(FitParamsPions[i][0])<<std::endl;
        //         std::cout<<std::get<2>(FitParamsPions[i][0])<<"*z + "<<std::get<3>(FitParamsPions[i][0])<<std::endl;
        //         std::cout<<"Photon2: "<<std::endl;
        //         std::cout<<std::get<0>(FitParamsPions[i][1])<<"*z + "<<std::get<1>(FitParamsPions[i][1])<<std::endl;
        //         std::cout<<std::get<2>(FitParamsPions[i][1])<<"*z + "<<std::get<3>(FitParamsPions[i][1])<<std::endl;
        //         std::cout<<"--------------------------------"<<std::endl;
        // }

        FindClosestApproach();
        TCanvas * Delta1 = new TCanvas("Delta1");
        Delta1->Divide(2,2,0.01,0.01);
        gStyle->SetOptStat(2222);
        TH1D * appdist1=new TH1D("Closest Approach", "ClosestApproach", 500,0,500);

        for(Int_t j = 0; j<nofEntries; j++) {
                //std::cout<<"ClosestApproach at: "<<std::get<0>(ClosestApproach[j]) <<"With distance: "<<std::get<1>(ClosestApproach[j])<<std::endl;
                appdist1->Fill(std::get<1>(ClosestApproach[j]));
        }
        appdist1->GetXaxis()->SetTitle("Distance of closest Approach [mm]");
        appdist1->GetYaxis()->SetTitle("Entries");


        CalculateDeviation();

        Delta1->cd(1);
        delx->Draw();
        Delta1->cd(2);
        dely->Draw();
        Delta1->cd(3);
        delz->Draw();

        Delta1->cd(4);
        appdist1->Draw();


        FitParamsPions.clear();
        ClosestApproach.clear();
}


//--------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::CalculateDeviation(){




        Double_t GunX=0.;
        Double_t GunY=0.;
        Double_t GunZ=-495.;

        for(Int_t i=0; i<nofEntries; i++) {

                Double_t ReconstX= ((std::get<0>(FitParamsPions[i][0])*std::get<0>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0])) +
                                    (std::get<0>(FitParamsPions[i][1])*std::get<0>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1])))/2;

                Double_t ReconstY= ((std::get<2>(FitParamsPions[i][0])*std::get<0>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]))+
                                    (std::get<2>(FitParamsPions[i][1])*std::get<0>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1])))/2;

                Double_t ReconstZ= std::get<0>(ClosestApproach[i]);



                Double_t DeltaX=GunX - ReconstX;
                Double_t DeltaY=GunY - ReconstY;
                Double_t DeltaZ=GunZ - ReconstZ;

                delx->Fill(DeltaX);
                dely->Fill(DeltaY);
                delz->Fill(DeltaZ);
        }

        delx->GetXaxis()->SetTitle("DeltaX[mm]");
        dely->GetXaxis()->SetTitle("DeltaY[mm]");
        delz->GetXaxis()->SetTitle("DeltaZ[mm]");

        delx->GetXaxis()->SetTitle("Entries");
        dely->GetYaxis()->SetTitle("Entries");
        delz->GetZaxis()->SetTitle("Entries");

        // Delta1->cd(1);
        // delx->Draw();
        //
        // Delta1->cd(2);
        // dely->Draw();
        //
        // Delta1->cd(3);
        // delz->Draw();





}

Double_t TROOTAnalysis::FindClosestApproach(){

        Fcn myfcnPi;

        myfcnPi.SetParamsPions(FitParamsPions);
        myfcnPi.SetApproachMode();

        MnUserParameters uparPi;

        double error_minimizer_parameters = 1e-4;

        uparPi.Add("z", 0., error_minimizer_parameters);

        MnMigrad migradPi(myfcnPi, uparPi, 2);

        for(Int_t i=0; i<nofEntries; i++) {

                uparPi.SetValue("z", 0.);

                myfcnPi.SetCurrentEvent(i);

                FunctionMinimum minPi = migradPi();

                //std::cout<<minPi<<std::endl;

                MnUserParameterState userParameterStatePi = minPi.UserState();

                //ClosestApproach.push_back(userParameterStatePi.Value("z"));

                Double_t X_=((std::get<0>(FitParamsPions[i][0])*userParameterStatePi.Value("z") + std::get<1>(FitParamsPions[i][0]))-
                             (std::get<0>(FitParamsPions[i][1])*userParameterStatePi.Value("z") + std::get<1>(FitParamsPions[i][1])));

                Double_t Y_=((std::get<2>(FitParamsPions[i][0])*userParameterStatePi.Value("z") + std::get<3>(FitParamsPions[i][0]))-
                             (std::get<2>(FitParamsPions[i][1])*userParameterStatePi.Value("z") + std::get<3>(FitParamsPions[i][1])));

                Double_t approachdist= TMath::Sqrt(X_*X_ + Y_*Y_);

                ClosestApproach.push_back(std::make_tuple(userParameterStatePi.Value("z"), approachdist));


        }

}


//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::CalcCOG(Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

        //variables for clustering
        Int_t xdir=1; // integers for direction
        Int_t ydir=0;
        Int_t buf;
        Int_t stepstodo=1; // how many steps to go before rotation
        Int_t stepsctr=0; // number of steps since last rotation

        Int_t currX=0; // starting the spiral on the histogram maximum
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

        for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) {    //loop over all simulated events
                std::cout<<"CenterEvent: "<<eventstodo<<std::endl;
                std::cout<<eventstodo-minevent<<std::endl;
                EcalTree->GetEntry(eventstodo); //grab evetn from tree
                Eges = Cevent->GapEnergy();
                Int_t cnh = Cevent->NHits();
                Double_t integSum=0;
                Double_t integral;
                Int_t i = showerCenters[eventstodo-minevent]-7;
                std::cout<<i<<std::endl;
                Int_t collected=0;
                if(i<0) {
                        i=0;
                        collected=4;
                }


                while(collected<7 && i<50) { //loop over layers until some layers after showerstart
                        if(i==49 && collected==0) {i=0; }
                        std::cout<<"i: "<<i<<std::endl;
                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i
                                if(Cevent->Hit(j)->Z()==i)
                                        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit()); //fill histogram with hits of layer i
                        }

                        integral=h1->Integral();

                        maxbin=h1->GetMaximumBin();
                        maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
                        h1->GetBinXYZ(maxbin, binx, biny, binz);

                        //Do spiral for clustering only if energy in Layer

                        if(integral!=0) {     //check if layer containes energy
                                collected++;
                                std::cout<<collected<<":"<<i<<std::endl;
                                xdir=1;        // integers for direction
                                ydir=0;
                                stepstodo=1;   // how many steps to go before rotation
                                stepsctr=0;   // number of steps since last rotation

                                currX=binx;    // starting the spiral on the histogram maximum
                                currY=biny;

                                esum=0;
                                curre=0;
                                esum+=maxdep;
                                h2->SetBinContent(currX, currY, maxdep);
                                auto tp=std::make_tuple(currX, currY, maxdep);
                                ClusteredHits.push_back(tp);
                                h1->SetBinContent(binx, biny, 0);

                                while(esum<integral*0.9 || currX<binx+4) {
                                        //do spiral until desired energyfraction is reached
                                        currX+=xdir;
                                        currY+=ydir;
                                        curre=h1->GetBinContent(currX, currY);

                                        if(curre>0) {     //save tile only if it containes energy
                                                tp=std::make_tuple(currX, currY, curre);
                                                ClusteredHits.push_back(tp);
                                                h2->SetBinContent(currX, currY, curre);
                                        }
                                        h1->SetBinContent(currX, currY, 0);
                                        esum+=curre;
                                        stepsctr++;
                                        if(stepsctr==stepstodo) {
                                                stepsctr=0;
                                                buf=xdir; //rotate 90 deg
                                                xdir= -ydir;
                                                ydir=buf;

                                                if(ydir==0) { //incremant steps at every iteration along x
                                                        stepstodo++;
                                                }
                                        }
                                }

                                Double_t e=h2->Integral(); //get energy contained in cluster

                                Eweight=e/Eges; //calculate weight

                                cgx=h2->GetMean(1);
                                xerr=h2->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h2->GetMean(2);
                                yerr=h2->GetMeanError(2);

                                if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                                if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                                cgz=i;

                                er1->Fill(i, xerr/nofEntries);
                                er2->Fill(i, yerr/nofEntries);

                                //h3->Fill(cgx, cgy, cgz);

                                auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
                                coglist.push_back(cg);

                        } // end of clustering

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::CalcCOGwithFit(Int_t minevent, Int_t maxevent){


        findShowercenter(minevent, maxevent);
        //CalcCOG(minevent, maxevent);                    //overloaded CalcCOG function for getting COGs after cluster start
        std::cout<<"center done"<<std::endl;
        //FitCOGs(minevent, maxevent, 5);
        std::cout<<"center fit done"<<std::endl;
        //PrintFitHists();
        COGCollection.clear(); // clear the cluster and COG vectors for reclustering
        coglist.clear();
        ClusteredHits.clear();


        //variables for clustering
        Int_t xdir=1; // integers for direction
        Int_t ydir=0;
        Int_t buf;
        Int_t stepstodo=1; // how many steps to go before rotation
        Int_t stepsctr=0; // number of steps since last rotation

        Int_t currX=0; // starting the spiral on the histogram maximum
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

        for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) {  //loop over all simulated events
                std::cout<<"ReClustering event: "<<eventstodo<<std::endl;
                EcalTree->GetEntry(eventstodo); //grab evetn from tree
                Eges = Cevent->GapEnergy();
                //std::cout<<Eges<<std::endl;
                Int_t cnh = Cevent->NHits();
                Double_t integSum=0;
                Double_t integral;


                for(Int_t i=0; i<50; i++) {
                        //  std::cout<<"layer: "<<i<<std::endl;
                        Int_t nofFills=0;             //loop over all layers in event
                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in event i
                                if(Cevent->Hit(j)->Z()==i) {
                                        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit()); //fill histogram with hits of layer i
                                        nofFills++;
                                }
                        }
                        //std::cout<<nofFills<<std::endl;

                        Double_t MoliereRaduis=48; //mm

                        Int_t nofLoops=0;

                        while(nofLoops<nofFills) {                   // check if maximum is within 1 MoliereRaduis of Fit

                                maxbin=h1->GetMaximumBin();
                                maxdep=h1->GetBinContent(maxbin);    // get coordinates of bin with maximum energy deposition
                                h1->GetBinXYZ(maxbin, binx, biny, binz);
                                Double_t Xpred=(std::get<0>(FitParamsGamma[eventstodo - minevent])*i+std::get<1>(FitParamsGamma[eventstodo - minevent]));
                                Double_t Ypred=(std::get<2>(FitParamsGamma[eventstodo - minevent])*i+std::get<3>(FitParamsGamma[eventstodo - minevent]));

                                //  std::cout<<"Maxdep: "<<binx<<":"<<biny<<"Pred: "<<Xpred<<":"<<Ypred<<std::endl;

                                Double_t d = TMath::Sqrt(    ((-495 + binx*10) - Xpred) * ((-495 + binx*10) - Xpred)  +   ((-495 + biny*10) - Ypred) * ((-495 + biny*10)-Ypred)   );

                                //  std::cout<<"d: "<<d<<std::endl;

                                if(d>MoliereRaduis) {

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

                        if(integral!=0) { //check if layer containes energy

                                xdir=1; // integers for direction
                                ydir=0;
                                stepstodo=1; // how many steps to go before rotation
                                stepsctr=0; // number of steps since last rotation

                                currX=binx; // starting the spiral on the histogram maximum
                                currY=biny;

                                esum=0;
                                curre=0;
                                esum+=maxdep;
                                h2->SetBinContent(currX, currY, maxdep);
                                auto tp=std::make_tuple(currX, currY, maxdep);
                                ClusteredHits.push_back(tp);
                                h1->SetBinContent(binx, biny, 0);

                                while(esum<integral*0.9 || currX<binx+9) {
                                        //do spiral until desired energyfraction is reached
                                        currX+=xdir;
                                        currY+=ydir;

                                        curre=h1->GetBinContent(currX, currY);

                                        if(curre>0) {     //save tile only if it containes energy
                                                tp=std::make_tuple(currX, currY, curre);
                                                ClusteredHits.push_back(tp);
                                                h2->SetBinContent(currX, currY, curre);
                                                h1->SetBinContent(currX, currY, 0);
                                        }


                                        esum+=curre;
                                        stepsctr++;
                                        if(stepsctr==stepstodo) {
                                                stepsctr=0;
                                                buf=xdir; //rotate 90 deg
                                                xdir= -ydir;
                                                ydir=buf;

                                                if(ydir==0) { //incremant steps at every iteration along x
                                                        stepstodo++;
                                                }
                                        }
                                }



                                Double_t e=h2->Integral(); //get energy contained in cluster

                                Eweight=e/Eges; //calculate weight

                                cgx=h2->GetMean(1);
                                xerr=h2->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h2->GetMean(2);
                                yerr=h2->GetMeanError(2);

                                if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                                if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                                cgz=i;

                                er1->Fill(i, xerr/nofEntries);
                                er2->Fill(i, yerr/nofEntries);



                                //std::cout<<"Center of gravity: "<<cgx<<" "<<cgy<<" "<<cgz<<std::endl;
                                h3->Fill(cgx, cgy, cgz);

                                auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
                                coglist.push_back(cg);

                        } // end of clustering


                        integSum += integral;



                        ClusteredHits.clear();

                        h1->Reset();

                        h2->Reset();

                        Eweight=0;




                }     // loop over layers

                COGCollection.push_back(coglist);
                coglist.clear();

        }             // loop over events
        FitParamsGamma.clear();
        //FitCOGs(minevent, maxevent, 5);
        PrintFitHists(minevent, maxevent);
        plotCOGs();
        //std::cout<<breakctr<<std::endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------





//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::FitCOGsPion( Int_t minevent, Int_t maxevent, Double_t tileLen, Bool_t isPion, Int_t photonNr){
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




        //transform from copynumber coordinates to geant4 coordinates
        std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > Transcoglist;
        std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > TransfomedCOGs;

        for(Int_t i = 0; i<nofEntries; i++) {                                 //transformation to geant4 coordinate system
                for(Int_t j = 0; j<COGCollection[i].size(); j++) {

                        auto tc = std::make_tuple(((-500+tileLen/2) + std::get<0>(COGCollection[i][j]) * tileLen),
                                                  ((-500+tileLen/2) + std::get<1>(COGCollection[i][j]) * tileLen),
                                                  -288.2 + std::get<2>(COGCollection[i][j]) * 11.8,
                                                  std::get<3>(COGCollection[i][j]),
                                                  std::get<4>(COGCollection[i][j]),
                                                  std::get<5>(COGCollection[i][j])  );
                        Transcoglist.push_back(tc);

                }
                TransfomedCOGs.push_back(Transcoglist);
                Transcoglist.clear();
        }

        // for(Int_t i = 0; i<nofEntries;i++){
        //   for(Int_t j =0;j<TransfomedCOGs[i].size();j++){
        //     std::cout<<std::get<0>(TransfomedCOGs[i][j])<<":"<<std::get<1>(TransfomedCOGs[i][j])<<":"<<std::get<2>(TransfomedCOGs[i][j])<<"-----"
        //     <<std::get<0>(COGCollection[i][j])<<":"<<std::get<1>(COGCollection[i][j])<<":"<<std::get<2>(COGCollection[i][j])<<std::endl;;
        //   }
        // }


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




        for(int events = 0; events < maxevent-minevent; events++) {

                upar.SetValue("mx", 0.);
                upar.SetValue("tx", 0.);
                upar.SetValue("my", 0.);
                upar.SetValue("ty", 0.);

                myfcn.SetCurrentEvent(events);

                std::cout<<"event: "<<events<<std::endl;
                FunctionMinimum min = migrad();

                std::cout<<min<<std::endl;


                MnUserCovariance cov = min.UserCovariance();

                //std::cout<<cov(0 ,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<cov(3,3)<<std::endl;


                MnUserParameterState userParameterState = min.UserState();




                if(isPion==false) {
                        auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                                                  userParameterState.Value("my"), userParameterState.Value("ty"));

                        FitParamsGamma.push_back(tp);
                }
                else if(isPion==true && photonNr==1) {
                        //std::cout<<"attempt adding params for photon 1"<<std::endl;

                        auto tp1 = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                                                   userParameterState.Value("my"), userParameterState.Value("ty"));
                        //std::cout<<"errga"<<std::endl;
                        FitParamsPions[events].push_back(tp1);
                        std::cout<<"adding params for photon 1"<<std::endl;
                }
                else if(isPion==true && photonNr==2) {
                        //std::cout<<"attempt adding params for photon 2"<<std::endl;
                        auto tp2 = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                                                   userParameterState.Value("my"), userParameterState.Value("ty"));

                        FitParamsPions[events].push_back(tp2);
                        std::cout<<"adding params for photon 2"<<std::endl;
                }

                Double_t covar[4][4];
                Double_t error[4];

                error[0]= userParameterState.Error("mx");
                error[1]= userParameterState.Error("tx");
                error[2]= userParameterState.Error("my");
                error[3]= userParameterState.Error("ty");

                //Fill correlation Histograms
                if(COGCollection[events].size() > 1) {


                        for(Int_t row=0; row<4; row++) {
                                for(Int_t col=0; col<4; col++) {
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

        }



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

        TransfomedCOGs.clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------



void TROOTAnalysis::PrintFitHists(Int_t minevent, Int_t maxevent){
        Int_t ent=maxevent-minevent;
        TCanvas * c1 = new TCanvas("Fit", "Fit", 1500,1300);

        c1->Divide(2,3,0.01, 0.01);

        //gStyle->SetOptFit(111111111111);


        for(Int_t i=0; i<ent; i++) {
                fitX->Fill(std::get<1>(FitParamsGamma[i]));
        }

        for(Int_t i=0; i<ent; i++) {
                fitY->Fill(std::get<3>(FitParamsGamma[i]));
        }

        for(Int_t i=0; i<ent; i++) {
                fitSX->Fill(std::get<0>(FitParamsGamma[i]));
        }

        for(Int_t i=0; i<ent; i++) {
                fitSY->Fill(std::get<2>(FitParamsGamma[i]));
        }

        c1->cd(1);
        fitX->GetXaxis()->SetTitle("X[mm]");
        fitX->GetYaxis()->SetTitle("Counts");
        gStyle->SetOptStat(1111);
        // Set stat options
        gStyle->SetStatY(0.9);
        // Set y-position (fraction of pad size)
        gStyle->SetStatX(0.9);
        // Set x-position (fraction of pad size)
        gStyle->SetStatW(0.25);
        // Set width of stat-box (fraction of pad size)
        gStyle->SetStatH(0.25);
        // Set height of stat-box (fraction of pad size)
        fitX->Draw();

        c1->cd(2);
        fitY->GetXaxis()->SetTitle("Y[mm]");
        fitY->GetYaxis()->SetTitle("Counts");
        fitY->Draw();

        c1->cd(3);
        fitSX->GetXaxis()->SetTitle("X Slope");
        fitSX->GetYaxis()->SetTitle("Counts");
        fitSX->Draw();

        c1->cd(4);
        fitSY->GetXaxis()->SetTitle("Y Slope");
        fitSY->GetYaxis()->SetTitle("Counts");
        fitSY->Draw();

        c1->cd(5);
        er1->Draw();

        c1->cd(6);
        er2->Draw();

        if(pathset) {

                std::string imgname= "Fit.png";
                std::string imgpath = savepath + '/' + imgname;
                TImage * img2 =TImage::Create();
                img2->FromPad(c1);
                img2->WriteImage(imgpath.c_str());
                delete img2;

        }
//  er1->Reset();
//  er2->Reset();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------


void TROOTAnalysis::CleanCOGs(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){



        TCanvas * c3 = new TCanvas("CleanCOGs", "CleanCOGs");
        TH3D * clCOG = new TH3D("COGs after cleaning", "COGs after celaning",20,40,60,20,40,60,50,0,50);

        Double_t MoliereRaduis=4.87; // in cm

        //CalcCOG(0,10, minevent-1, maxevent-1);
        //FitCOGs(minevent-1, maxevent-1, 5);

        COGCollection.clear();
        coglist.clear();
        ClusteredHits.clear();

        //PrintFitHists();
        //CalcCOG(minlayer-1,maxlayer-1,minevent-1, maxevent-1 );

        Int_t erasecounter=0;
        Int_t totalcounter=0;
        for(int i=0; i<nofEntries; i++) {                               //Delete the COGs, outside of one MoliereRadius
                for(int z=0; z<COGCollection[i].size(); z++) {

                        Double_t distx=(std::get<0>(FitParamsGamma[i])*z+std::get<1>(FitParamsGamma[i]))-std::get<0>(COGCollection[i][z]);
                        Double_t disty=(std::get<2>(FitParamsGamma[i])*z+std::get<3>(FitParamsGamma[i]))-std::get<1>(COGCollection[i][z]);
                        totalcounter++;
                        if(TMath::Sqrt(distx*distx)>MoliereRaduis || TMath::Sqrt(disty*disty)>MoliereRaduis) {
                                COGCollection[i].erase(COGCollection[i].begin()+z);
                                //std::cout<<"X distance: "<<distx<<" Y distance: "<<disty<<std::endl;
                                erasecounter++;
                        }

                }

        }

        for(int i=0; i<nofEntries; i++) {                               //Delete the COGs, outside of one MoliereRadius
                for(int z=0; z<COGCollection[i].size(); z++) {
                        clCOG->Fill(std::get<0>(COGCollection[i][z]), std::get<1>(COGCollection[i][z]), std::get<2>(COGCollection[i][z]));
                }
        }

        clCOG->GetXaxis()->SetTitle("X");
        clCOG->GetYaxis()->SetTitle("Y");
        clCOG->GetZaxis()->SetTitle("Z");
        clCOG->SetMarkerStyle(4);
        c3->cd();
        clCOG->Draw("");

        FitParamsGamma.clear();


        //FitCOGs(minevent, maxevent, 5);          //refit the leftover COGs
        //PrintFitHists(minevent, maxevent);
        //plotCOGs();

        std::cout<<"Discarded COGs: "<<erasecounter<<" Total COGs:"<<totalcounter<<std::endl;




}

//////////////////////////////Old Section///////////////////////////////////

void TROOTAnalysis::CalcCOG(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){            //calculate vector of (X,Y) tuples containing layerwise center of gravity

        //variables for clustering
        Int_t xdir=1; // integers for direction
        Int_t ydir=0;
        Int_t buf;
        Int_t stepstodo=1; // how many steps to go before rotation
        Int_t stepsctr=0; // number of steps since last rotation

        Int_t currX=0; // starting the spiral on the histogram maximum
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

        for(Int_t eventstodo=minevent; eventstodo<maxevent; eventstodo++) {    //loop over all simulated events

                EcalTree->GetEntry(eventstodo);       //grab evetn from tree
                Eges = Cevent->GapEnergy();
                Int_t cnh = Cevent->NHits();
                Double_t integral;
                Bool_t foundstart=false;
                for(Int_t i=minlayer; i<maxlayer; i++) { //loop over all layers in event
                        for(Int_t j =0; j<cnh; j++) { //loop over all hits in laver i
                                if(Cevent->Hit(j)->Z()==i)
                                        h1->Fill(Cevent->Hit(j)->X(), Cevent->Hit(j)->Y(), Cevent->Hit(j)->EnergyDeposit());
                                //fill histogram with hits of layer i
                        }

                        integral=h1->Integral();

                        maxbin=h1->GetMaximumBin();
                        maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
                        h1->GetBinXYZ(maxbin, binx, biny, binz);

                        //Do spiral for clustering only if energy in Layer

                        if(integral!=0) {     //check if layer containes energy

                                xdir=1;        // integers for direction
                                ydir=0;
                                stepstodo=1;   // how many steps to go before rotation
                                stepsctr=0;   // number of steps since last rotation

                                currX=binx;    // starting the spiral on the histogram maximum
                                currY=biny;

                                esum=0;
                                curre=0;
                                esum+=maxdep;
                                h2->SetBinContent(currX, currY, maxdep);
                                auto tp=std::make_tuple(currX, currY, maxdep);
                                ClusteredHits.push_back(tp);
                                h1->SetBinContent(binx, biny, 0);

                                while(esum<integral*0.9 || currX<binx+4) {
                                        //do spiral until desired energyfraction is reached
                                        currX+=xdir;
                                        currY+=ydir;
                                        curre=h1->GetBinContent(currX, currY);

                                        if(curre>0) {     //save tile only if it containes energy
                                                tp=std::make_tuple(currX, currY, curre);
                                                ClusteredHits.push_back(tp);
                                                h2->SetBinContent(currX, currY, curre);
                                        }
                                        h1->SetBinContent(currX, currY, 0);
                                        esum+=curre;
                                        stepsctr++;
                                        if(stepsctr==stepstodo) {
                                                stepsctr=0;
                                                buf=xdir; //rotate 90 deg
                                                xdir= -ydir;
                                                ydir=buf;

                                                if(ydir==0) { //incremant steps at every iteration along x
                                                        stepstodo++;
                                                }
                                        }
                                }

                                Double_t e=h2->Integral(); //get energy contained in cluster

                                Eweight=e/Eges; //calculate weight

                                cgx=h2->GetMean(1);
                                xerr=h2->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h2->GetMean(2);
                                yerr=h2->GetMeanError(2);

                                if(yerr<0.00001) {yerr=1/TMath::Sqrt(12); }
                                if(xerr<0.00001) {xerr=1/TMath::Sqrt(12); }

                                cgz=i;

                                er1->Fill(i, xerr/nofEntries);
                                er2->Fill(i, yerr/nofEntries);

                                h3->Fill(cgx, cgy, cgz);

                                auto cg=std::make_tuple(cgx,cgy,cgz, xerr, yerr, Eweight);
                                coglist.push_back(cg);

                        } // end of clustering

                        ClusteredHits.clear();

                        h1->Reset();

                        h2->Reset();

                }
                COGCollection.push_back(coglist);
                coglist.clear();

        }
        plotCOGs();
}


void TROOTAnalysis::FitCOGs( Int_t minevent, Int_t maxevent, Double_t tileLen){
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




        //transform from copynumber coordinates to geant4 coordinates
        std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > Transcoglist;
        std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > TransfomedCOGs;

        for(Int_t i = 0; i<nofEntries; i++) {                                 //transformation to geant4 coordinate system
                for(Int_t j = 0; j<COGCollection[i].size(); j++) {

                        auto tc = std::make_tuple(((-500+tileLen/2) + std::get<0>(COGCollection[i][j]) * tileLen),
                                                  ((-500+tileLen/2) + std::get<1>(COGCollection[i][j]) * tileLen),
                                                  -288.2 + std::get<2>(COGCollection[i][j]) * 11.8,
                                                  std::get<3>(COGCollection[i][j]),
                                                  std::get<4>(COGCollection[i][j]),
                                                  std::get<5>(COGCollection[i][j])  );
                        Transcoglist.push_back(tc);

                }
                TransfomedCOGs.push_back(Transcoglist);
                Transcoglist.clear();
        }

        // for(Int_t i = 0; i<nofEntries;i++){
        //   for(Int_t j =0;j<TransfomedCOGs[i].size();j++){
        //     std::cout<<std::get<0>(TransfomedCOGs[i][j])<<":"<<std::get<1>(TransfomedCOGs[i][j])<<":"<<std::get<2>(TransfomedCOGs[i][j])<<"-----"
        //     <<std::get<0>(COGCollection[i][j])<<":"<<std::get<1>(COGCollection[i][j])<<":"<<std::get<2>(COGCollection[i][j])<<std::endl;;
        //   }
        // }


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



        for(int events = 0; events < maxevent-minevent; events++) {

                upar.SetValue("mx", 0.);
                upar.SetValue("tx", 0.);
                upar.SetValue("my", 0.);
                upar.SetValue("ty", 0.);

                myfcn.SetCurrentEvent(events);

                std::cout<<"event: "<<events<<std::endl;
                FunctionMinimum min = migrad();

                std::cout<<min<<std::endl;


                MnUserCovariance cov = min.UserCovariance();

                //std::cout<<cov(0 ,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<cov(3,3)<<std::endl;


                MnUserParameterState userParameterState = min.UserState();

                auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                                          userParameterState.Value("my"), userParameterState.Value("ty"));
                FitParamsGamma.push_back(tp);



                Double_t covar[4][4];
                Double_t error[4];

                error[0]= userParameterState.Error("mx");
                error[1]= userParameterState.Error("tx");
                error[2]= userParameterState.Error("my");
                error[3]= userParameterState.Error("ty");

                //Fill correlation Histograms
                if(COGCollection[events].size() > 1) {


                        for(Int_t row=0; row<4; row++) {
                                for(Int_t col=0; col<4; col++) {
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
