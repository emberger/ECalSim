#include "TROOTAnalysis.hh"




TROOTAnalysis::TROOTAnalysis(TChain* ch){
        this->EcalTree=ch->GetTree();
        std::cout << "assigned Tree" << '\n';
        nofEntries=EcalTree->GetEntries();
        std::cout<<nofEntries<<std::endl;
        FitParamsPions.resize(nofEntries);


        Cevent = new B4ROOTEvent();
        EcalTree->SetBranchAddress("EventBranch", &Cevent);


        for(Int_t i=0; i<nofEntries; i++) {
                EcalTree->GetEntry(i);
                if(i==0) {
                        GapThickness=Cevent->GapThickness();
                        AbsoThickness=Cevent->AbsoThickness();
                        tiledimX=Cevent->TilesizeX();
                        tiledimY=Cevent->TilesizeY();
                        calsizeXY=Cevent->calsizeXY();
                        nofLayers=Cevent->NumberOfLayers();
                }
                Int_t hitnr = Cevent->NHits();
                Double_t tmpE1=0;
                Double_t tmpE2=0;
                for(Int_t j=0; j<hitnr; j++) {

                        if(Cevent->Hit(j)->PhotNr()==1) {
                                tmpE1+=Cevent->Hit(j)->EnergyDeposit();
                        }
                        else if(Cevent->Hit(j)->PhotNr()==2) {
                                tmpE2+=Cevent->Hit(j)->EnergyDeposit();
                        }

                }
                EnergyPhoton1.push_back(tmpE1/0.50);
                EnergyPhoton2.push_back(tmpE2/0.50);

        }

        histsizeX=calsizeXY/tiledimX;
        histsizeY=calsizeXY/tiledimY;
        histsizeZ=nofLayers;

        h = new TH3D("ECalEvent","ECalEvent",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);


        hA = new TH3D("ECalEvent1","ECalEvent1",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
        hB = new TH3D("ECalEvent2","ECalEvent2",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);
        h1 = new TH2D("h1", "h1", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);
        h2 = new TH2D("h2", "h2", histsizeX+1,-0.5,histsizeX+0.5,histsizeY+1,-0.5,histsizeY+0.5);
        h3 = new TH3D("h3", "h3",histsizeX,0,histsizeX,histsizeY,0,histsizeY,histsizeZ,0,histsizeZ);

        std::cout<<"GT: "<<GapThickness<<" AT: "<<AbsoThickness<<" TDX: "<<tiledimX<<" TDY: "<<tiledimY<<" CXY: "<<calsizeXY <<" NL: "<<nofLayers<<std::endl;
        std::cout<<" HX: "<<histsizeX<<" HY: "<<histsizeY<<" HZ: "<<histsizeZ<<std::endl;
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

        std::string Path="/home/iwsatlas1/emberger/PionPlots/Event";
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
        rms->Divide(2,2,0.01,0.01);



        // Double_t y[10]={0.5298,0.3507,0.3624,0.2524,0.222,0.1861,0.228,0.1751,0.1425,0.1481}; // 0.0 slope
        // Double_t x[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        //
        Double_t x1[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        Double_t y1[10]={64.25/102,55.13/102,37.8/102,37.7/102,36.5/102,38.6/102,23.0/102,25.7/102,16.3/102,30.5/102};            //0.2slope 10cm

        Double_t x2[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        Double_t y2[10]={76.2/107,53.5/107,45.8/107,63.6/107,56.1/107,42.1/107,34.5/107,27.68/107,35.1/107,28.33/107};            //0.4slope 10cm

        Double_t x3[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        Double_t y3[10]={100.4/203,62.34/203,52.65/203,54.66/203,48.39/203,40.51/203,48.15/203,35.07/203,36.26/203,39.16/203};            //0.2slope 20cm

        Double_t x4[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
        Double_t y4[10]={105.8/215,80.74/215,67.6/215,66.6/215,58.37/215,42.73/215,51.5/215,45.47/215,47.5/215,42.97/215};            //0.4slope 20cm

        TF1 * rms_x = new TF1("rms", "sqrt(([0] * [0] / x) + ([2]*[2])  + ([1]*[1]/(x*x)))");

        TGraph *gRMS1=new TGraph(10, x1, y1 );
        TGraph *gRMS2=new TGraph(10, x2, y2 );
        TGraph *gRMS3=new TGraph(10, x3, y3 );
        TGraph *gRMS4=new TGraph(10, x4, y4 );

        gRMS1->Fit(rms_x);
        gStyle->SetOptFit(11111111);
        gRMS1->SetTitle("Distance: 10cm , slope 0.2");
        gRMS1->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS1->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS1->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(1);
        gRMS1->Draw("A*");

        gRMS2->Fit(rms_x);
        gRMS2->SetTitle("Distance: 10cm , slope 0.4");
        gRMS2->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS2->GetYaxis()->SetTitle("#frac{RMS_x}{100mm}");
        gRMS2->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(2);
        gRMS2->Draw("A*");

        gRMS3->Fit(rms_x);
        gRMS3->SetTitle("Distance: 20cm , slope 0.2");
        gRMS3->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS3->GetYaxis()->SetTitle("#frac{RMS_x}{200mm}");
        gRMS3->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(3);
        gRMS3->Draw("A*");

        gRMS4->Fit(rms_x);
        gRMS4->SetTitle("Distance: 20cm , slope 0.4");
        gRMS4->GetXaxis()->SetTitle("Energy[GeV]");
        gRMS4->GetYaxis()->SetTitle("#frac{RMS_x}{200mm}");
        gRMS4->GetYaxis()->SetTitleOffset(1.3);
        rms->cd(4);
        gRMS4->Draw("A*");

        rms->Print("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/rmsx_over_E_analysis/Res.png");
        rms->Print("/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/rmsx_over_E_analysis/Res.C");


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
                Double_t proX=std::get<0>(FitParamsGamma[i])*(-275-distance)+std::get<1>(FitParamsGamma[i]);
                Double_t proY=std::get<2>(FitParamsGamma[i])*(-275-distance)+std::get<3>(FitParamsGamma[i]);

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
                        std::cout<<"doing event"<<eventstodo<<", photon "<<photonNR<<std::endl;
                        maxbin=h1->GetMaximumBin();
                        maxdep=h1->GetBinContent(maxbin); // get coordinates of bin with maximum energy deposition
                        h1->GetBinXYZ(maxbin, binx, biny, binz);

                        //Do spiral for clustering only if energy in Layer

                        if(integral!=0) { //check if layer containes energy

                                // xdir=1;  // integers for direction
                                // ydir=0;
                                // stepstodo=1; // how many steps to go before rotation
                                // stepsctr=0; // number of steps since last rotation
                                //
                                // currX=binx; // starting the spiral on the histogram maximum
                                // currY=biny;
                                //
                                // esum=0;
                                // curre=0;
                                // esum+=maxdep;
                                // h2->SetBinContent(currX, currY, maxdep);
                                // auto tp=std::make_tuple(currX, currY, maxdep);
                                // ClusteredHits.push_back(tp);
                                // h1->SetBinContent(binx, biny, 0);

                                // while(esum<integral*0.9 || currX<binx+4) {
                                //         //do spiral until desired energyfraction is reached
                                //         currX+=xdir;
                                //         currY+=ydir;
                                //         curre=h1->GetBinContent(currX, currY);
                                //
                                //         if(curre>0) { //save tile only if it containes energy
                                //                 tp=std::make_tuple(currX, currY, curre);
                                //                 ClusteredHits.push_back(tp);
                                //                 h2->SetBinContent(currX, currY, curre);
                                //         }
                                //         h1->SetBinContent(currX, currY, 0);
                                //         esum+=curre;
                                //         stepsctr++;
                                //         if(stepsctr==stepstodo) {
                                //                 stepsctr=0;
                                //                 buf=xdir; //rotate 90 deg
                                //                 xdir= -ydir;
                                //                 ydir=buf;
                                //
                                //                 if(ydir==0) { //incremant steps at every iteration along x
                                //                         stepstodo++;
                                //                 }
                                //         }
                                // }

                                Double_t e=h1->Integral(); //get energy contained in cluster

                                if(photonNR==1) {
                                        Eweight=e/Cevent->EnergyPhoton1(); //calculate weight
                                }
                                else if(photonNR==2) {
                                        Eweight=e/Cevent->EnergyPhoton2();
                                }

                                cgx=h1->GetMean(1);
                                xerr=h1->GetMeanError(1);
                                //extract center of gravity and error
                                cgy=h1->GetMean(2);
                                yerr=h1->GetMeanError(2);

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

                }         //end of event

                if(photonNR==1) {
                        showerCOGPhoton1.push_back(std::make_tuple(h3->GetMean(1), h3->GetMean(2), h3->GetMean(3)));
                        std::cout<<"111"<<std::endl;
                }

                else if(photonNR==2) {
                        showerCOGPhoton2.push_back(std::make_tuple(-1*(h3->GetMean(1)), -1*(h3->GetMean(2)), h3->GetMean(3)));
                        std::cout<<"222"<<std::endl;
                }

                h3->Reset();
                COGCollection.push_back(coglist);
                coglist.clear();

        }
        //plotCOGs();
}


//-------------------------------------------------------------------------------------------------------------------------------------
void TROOTAnalysis::FitCOGsPion( Int_t minevent, Int_t maxevent, Bool_t isPion, Int_t photonNr){
        //histograms for correlation
        // TCanvas * corr1 = new TCanvas("Correlations");
        // corr1->Divide(3,2,0.01,0.01);
        //
        // TH1D * co1 = new TH1D("MxTx correllation", "MxTx correllation",100, -1,1 );
        // co1->GetXaxis()->SetTitle("corellation of X slope and X intercept");
        // TH1D * co2 = new TH1D("MxMy correllation", "MxMy correllation",100, -1,1 );
        // co2->GetXaxis()->SetTitle("corellation of X slope and Y slope");
        // TH1D * co3 = new TH1D("MxTy correllation", "MxTy correllation",100, -1,1 );
        // co3->GetXaxis()->SetTitle("correlation of X slope and Y intercept");
        // TH1D * co4 = new TH1D("TxMy correllation", "TxMy correllation",100, -1,1 );
        // co4->GetXaxis()->SetTitle("correlation of X intercept and Y slope");
        // TH1D * co5 = new TH1D("TxTy correllation", "TxTy correllation",100, -1,1 );
        // co5->GetXaxis()->SetTitle("correlation of X intercapt an Y intercept");
        // TH1D * co6 = new TH1D("MyTy correllation", "MyTy correllation",100, -1,1 );
        // co6->GetXaxis()->SetTitle("correlation of Y slope and Y intercept");




        //transform from copynumber coordinates to geant4 coordinates
        std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > Transcoglist;
        std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > TransfomedCOGs;


//----------------------------------------------------------------1.8mm abso-----------------------------------------------------------------------------------------


        // for(Int_t i = 0; i<nofEntries; i++) {                                 //transformation to geant4 coordinate system
        //         for(Int_t j = 0; j<COGCollection[i].size(); j++) {
        //
        //                 auto tc = std::make_tuple(((-500+tileLen/2) + std::get<0>(COGCollection[i][j]) * tileLen),
        //                                           ((-500+tileLen/2) + std::get<1>(COGCollection[i][j]) * tileLen),
        //                                           -288.2 + std::get<2>(COGCollection[i][j]) * 11.8,
        //                                           std::get<3>(COGCollection[i][j]),
        //                                           std::get<4>(COGCollection[i][j]),
        //                                           std::get<5>(COGCollection[i][j])  );
        //                 Transcoglist.push_back(tc);
        //
        //         }
        //         TransfomedCOGs.push_back(Transcoglist);
        //         Transcoglist.clear();
        // }

//----------------------------------------------------------------1mm Abso-----------------------------------------------------------------------------------------
        Double_t offsetX=((-1*calsizeXY/2)+tiledimX/2);
        Double_t offsetY=((-1*calsizeXY/2)+tiledimY/2);
        Double_t offsetZ=((((AbsoThickness+GapThickness)*(-1)*nofLayers)/2)+(AbsoThickness+GapThickness/2));
        std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

        for(Int_t i = 0; i<nofEntries; i++) {                                         //transformation to geant4 coordinate system
                for(Int_t j = 0; j<COGCollection[i].size(); j++) {

                        auto tc = std::make_tuple((offsetX + std::get<0>(COGCollection[i][j]) * tiledimX),                // X
                                                  (offsetY + std::get<1>(COGCollection[i][j]) * tiledimY),                 // Y
                                                  (offsetZ + std::get<2>(COGCollection[i][j]) * (AbsoThickness+GapThickness)),                               // Z
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
        myfcn.SetMode("photon");
        myfcn.SetCOGs(TransfomedCOGs);
        //myfcn.PrintCOGs();
        MnUserParameters upar;
        double error_minimizer_parameters = 1e-4;

        upar.Add("mx", 0., error_minimizer_parameters);
        upar.Add("tx", 1., error_minimizer_parameters);
        upar.Add("my", 0., error_minimizer_parameters);
        upar.Add("ty", 1., error_minimizer_parameters);

        //cout << upar << endl;

        MnMigrad migrad(myfcn, upar, 2);




        for(int events = 0; events < nofEntries; events++) {

                upar.SetValue("mx", 0.);
                upar.SetValue("tx", 0.);
                upar.SetValue("my", 0.);
                upar.SetValue("ty", 0.);

                myfcn.SetCurrentEvent(events);

                //std::cout<<"event: "<<events<<std::endl;
                FunctionMinimum min = migrad();

                //std::cout<<min<<std::endl;


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

                //         Double_t covar[4][4];
                //         Double_t error[4];
                //
                //         error[0]= userParameterState.Error("mx");
                //         error[1]= userParameterState.Error("tx");
                //         error[2]= userParameterState.Error("my");
                //         error[3]= userParameterState.Error("ty");
                //
                //         //Fill correlation Histograms
                //         if(COGCollection[events].size() > 1) {
                //
                //
                //                 for(Int_t row=0; row<4; row++) {
                //                         for(Int_t col=0; col<4; col++) {
                //                                 covar[row][col]=cov(row,col);
                //
                //                         }
                //                 }
                //
                //                 co1->Fill(covar[0][1]/(error[0]*error[1]));
                //                 co2->Fill(covar[0][2]/(error[0]*error[2]));
                //                 co3->Fill(covar[0][3]/(error[0]*error[3]));
                //                 co4->Fill(covar[1][2]/(error[1]*error[2]));
                //                 co5->Fill(covar[1][3]/(error[1]*error[3]));
                //                 co6->Fill(covar[2][3]/(error[2]*error[3]));
                //         }
                //
        }
        //
        //
        //
        // corr1->cd(1);
        // co1->Draw("");
        //
        // corr1->cd(2);
        // co2->Draw("");
        //
        // corr1->cd(3);
        // co3->Draw("");
        //
        // corr1->cd(4);
        // co4->Draw("");
        //
        // corr1->cd(5);
        // co5->Draw("");
        //
        // corr1->cd(6);
        // co6->Draw("");

        TransfomedCOGs.clear();
}

//--------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::AnalyzePions(Int_t minlayer, Int_t maxlayer, Int_t minevent, Int_t maxevent){

        Bool_t pion=true;



        CalcCOGPion(minlayer, maxlayer, minevent, maxevent, 1);
        //plotCOGs();
        FitCOGsPion(minevent, maxevent, pion, 1);

        COGCollection.clear();
        coglist.clear();

        CalcCOGPion(minlayer, maxlayer, minevent, maxevent, 2);

        FitCOGsPion(minevent, maxevent, pion, 2);

        // TCanvas * lsk = new TCanvas("lsk");
        //
        //
        // TH3D * cogph1 = new TH3D("cogph1", "cogph1", 1000,-200,200,1000,-200,200,100,0,50);
        // TH3D * cogph2 = new TH3D("cogph2", "cogph2", 1000,-200,200,1000,-200,200,100,0,50);
        //
        // for(Int_t i =0; i<nofEntries; i++) {
        //         cogph1->Fill(std::get<0>(showerCOGPhoton1[i]),std::get<1>(showerCOGPhoton1[i]),std::get<2>(showerCOGPhoton1[i]));
        //         cogph2->Fill(std::get<0>(showerCOGPhoton2[i]),std::get<1>(showerCOGPhoton2[i]),std::get<2>(showerCOGPhoton2[i]));
        //         //         std::cout<<"Event: "<<i<<std::endl;
        //         //         std::cout<<std::get<0>(showerCOGPhoton1[i])<<":"<<std::get<1>(showerCOGPhoton1[i])<<":"<<std::get<2>(showerCOGPhoton1[i])<<std::endl;
        //         //         std::cout<<std::get<0>(showerCOGPhoton2[i])<<":"<<std::get<1>(showerCOGPhoton2[i])<<":"<<std::get<2>(showerCOGPhoton2[i])<<std::endl;
        // }
        //
        //
        //
        // cogph1->SetMarkerColor(kRed);
        // cogph1->SetMarkerStyle(20);
        // cogph1->SetMarkerSize(0.3);
        // cogph1->Draw();
        //
        // cogph2->SetMarkerColor(kBlue);
        // cogph2->SetMarkerStyle(20);
        // cogph2->SetMarkerSize(0.3);
        //
        // cogph2->Draw("same");

        FindClosestApproach();

        CalculateDeviation();

        GetInvariantMass();
        std::cout<<"PlottingChimap"<<std::endl;
        //PlotChimap(3);
        std::cout<<"DoneChimap"<<std::endl;


        //MinimizeClosestApproach();
        PionLocator();

        FitParamsPions.clear();
        ClosestApproach.clear();
}


//--------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::CalculateDeviation(){




        Double_t GunX=0.;
        Double_t GunY=0.;
        Double_t GunZ=-275.;

        for(Int_t i=0; i<nofEntries; i++) {







                Double_t DeltaX=GunX - std::get<0>(ClosestApproach[i]);;
                Double_t DeltaY=GunY - std::get<1>(ClosestApproach[i]);;
                Double_t DeltaZ=GunZ - std::get<2>(ClosestApproach[i]);;

                delx->Fill(DeltaX);
                dely->Fill(DeltaY);
                delz->Fill(DeltaZ);
        }

        TCanvas * Delta1 = new TCanvas("Delta1", "difference of gun positoin to ClosestApproach");
        Delta1->Divide(2,2,0.01,0.01);
        gStyle->SetOptStat(2222);


        delx->GetXaxis()->SetTitle("DeltaX[mm]");
        dely->GetXaxis()->SetTitle("DeltaY[mm]");
        delz->GetXaxis()->SetTitle("DeltaZ[mm]");

        delx->GetYaxis()->SetTitle("Entries");
        dely->GetYaxis()->SetTitle("Entries");
        delz->GetYaxis()->SetTitle("Entries");

        Delta1->cd(1);
        delx->Draw();
        Delta1->cd(2);
        dely->Draw();
        Delta1->cd(3);
        delz->Draw();

        for(Int_t j = 0; j<nofEntries; j++) {
                //std::cout<<"ClosestApproach at: "<<std::get<0>(ClosestApproach[j]) <<"With distance: "<<std::get<1>(ClosestApproach[j])<<std::endl;
                appdist1->Fill(std::get<3>(ClosestApproach[j]));
        }
        appdist1->GetXaxis()->SetTitle("Distance of closest Approach [mm]");
        appdist1->GetYaxis()->SetTitle("Entries");
        Delta1->cd(4);
        appdist1->Draw();



}

Double_t TROOTAnalysis::FindClosestApproach(){

        Fcn myfcnPi;

        myfcnPi.SetParamsPions(FitParamsPions);
        myfcnPi.SetMode("pion");

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


                Double_t mx1=(std::get<0>(FitParamsPions[i][0])*userParameterStatePi.Value("z") + std::get<1>(FitParamsPions[i][0]));
                Double_t mx2=(std::get<0>(FitParamsPions[i][1])*userParameterStatePi.Value("z") + std::get<1>(FitParamsPions[i][1]));

                Double_t my1=(std::get<2>(FitParamsPions[i][0])*userParameterStatePi.Value("z") + std::get<3>(FitParamsPions[i][0]));
                Double_t my2=(std::get<2>(FitParamsPions[i][1])*userParameterStatePi.Value("z") + std::get<3>(FitParamsPions[i][1]));

                Double_t ReconstX= (mx1 + mx2)/2;

                Double_t ReconstY= (my1 + my2)/2;

                Double_t X_=(mx1-mx2);

                Double_t Y_=(my1-my2);

                Double_t approachdist= TMath::Sqrt(X_*X_ + Y_*Y_);

                ClosestApproach.push_back(std::make_tuple(ReconstX, ReconstY, userParameterStatePi.Value("z"), approachdist));


        }

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void TROOTAnalysis::GetInvariantMass(){


        TVector3 V_ph1_0;                //vectors to construct direction vectors
        TVector3 V_ph2_0;

        TVector3 V_ph1_1;
        TVector3 V_ph2_1;


        TVector3 Ph1_direction;         //reconstructed direction unitvectors
        TVector3 Ph2_direction;

        TVector3 Ph1_dir_true;
        TVector3 Ph2_dir_true;

        TLorentzVector L_ph1;
        TLorentzVector L_ph2;
        TLorentzVector L_pi;

        TLorentzVector L_ph1_reco;
        TLorentzVector L_ph2_reco;
        TLorentzVector L_pi_reco;

        TLorentzVector L_ph1_true;
        TLorentzVector L_ph2_true;
        TLorentzVector L_pi_true;

        TCanvas * IM1 = new TCanvas("InvariantMassC");
        IM1->Divide(3,2,0.01,0.01);
        gStyle->SetOptStat(2222);


// Histograms for difference between true and reconstructed energy

        TCanvas * DirDelta1 = new TCanvas("DeltasDirectionC", "Difference of true and reconstructed photon direction");
        DirDelta1->Divide(2,3,0.01,0.01);

        TH1D * delX1 = new TH1D("DeltaX_Ph1", "DeltaXD_Ph1", 1000,-2,2);
        delX1->GetXaxis()->SetTitle("DeltaX");
        delX1->GetYaxis()->SetTitle("Entries");

        TH1D * delY1 = new TH1D("DeltaYD_Ph1", "DeltaYD_Ph1", 1000,-2,2);
        delY1->GetXaxis()->SetTitle("DeltaY");
        delY1->GetYaxis()->SetTitle("Entries");

        TH1D * delZ1 = new TH1D("DeltaZD_Ph1", "DeltaZD_Ph1", 1000,-2,2);
        delZ1->GetXaxis()->SetTitle("DeltaZ");
        delZ1->GetYaxis()->SetTitle("Entries");

        TH1D * delX2 = new TH1D("DeltaX_Ph2", "DeltaXD_Ph2", 1000,-2,2);
        delX2->GetXaxis()->SetTitle("DeltaX");
        delX2->GetYaxis()->SetTitle("Entries");

        TH1D * delY2 = new TH1D("DeltaYD_Ph2", "DeltaYD_Ph2", 1000,-2,2);
        delY2->GetXaxis()->SetTitle("DeltaY");
        delY2->GetYaxis()->SetTitle("Entries");

        TH1D * delZ2 = new TH1D("DeltaZD_Ph2", "DeltaZD_Ph2", 1000,-2,2);
        delZ2->GetXaxis()->SetTitle("DeltaZ");
        delZ2->GetYaxis()->SetTitle("Entries");


        // histograms for different invariaont mass calculations


        TH1D * InvMassReco = new TH1D("InvariantMassH1", "true Energy reconstructed direction",500, 0,1000);
        InvMassReco->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassReco->GetYaxis()->SetTitle("Entries");

        TH1D * InvMassSim = new TH1D("InvariantMassH2", "true Direction reconstructed Energy", 500,0,1000);
        InvMassSim->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassSim->GetYaxis()->SetTitle("Entries");

        TH1D * InvMassAllreco = new TH1D("InvariantMassH3", "reconstructed Direction reconstructed Energy", 500,0,1000);
        InvMassAllreco->GetXaxis()->SetTitle("InvariantMass[MeV]");
        InvMassAllreco->GetYaxis()->SetTitle("Entries");

        //histograms for Angle_vs_InvariantMass plots

        TH2D * Angle_vs_massReco = new TH2D("Angle_vs_massRecoH", "Angle_vs_InvariantMass", 1000,0,1000,500,0,200);
        Angle_vs_massReco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massReco->GetYaxis()->SetTitle("Angle[deg]");

        TH2D * Angle_vs_massSim = new TH2D("Angle_vs_massSimH", "Angle_vs_InvariantMass", 1000,0,1000,500,0,200);
        Angle_vs_massSim->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massSim->GetYaxis()->SetTitle("Angle[deg]");

        TH2D * Angle_vs_massAllreco = new TH2D("Angle_vs_massAllrecoH", "Angle_vs_InvariantMass", 1000,0,1000,500,0,200);
        Angle_vs_massAllreco->GetXaxis()->SetTitle("Invariant Mass [MeV]");
        Angle_vs_massAllreco->GetYaxis()->SetTitle("Angle[deg]");



        for(Int_t i=0; i<nofEntries; i++) {

                //calculate reconstructed direction vector Photon 1
                V_ph1_0.SetX(std::get<0>(FitParamsPions[i][0]) * 0 + std::get<1>(FitParamsPions[i][0]));
                V_ph1_0.SetY(std::get<2>(FitParamsPions[i][0]) * 0 + std::get<3>(FitParamsPions[i][0]));
                V_ph1_0.SetZ(0);

                V_ph1_1.SetX(std::get<0>(FitParamsPions[i][0]) * 1 + std::get<1>(FitParamsPions[i][0]));
                V_ph1_1.SetY(std::get<2>(FitParamsPions[i][0]) * 1 + std::get<3>(FitParamsPions[i][0]));
                V_ph1_1.SetZ(1);

                TVector3 V_ph1=V_ph1_1 - V_ph1_0;
                TVector3 Ph1_direction = V_ph1.Unit();

                //calculate reconstructed direction vector Photon 2
                V_ph2_0.SetX(std::get<0>(FitParamsPions[i][1]) * 0 + std::get<1>(FitParamsPions[i][1]));
                V_ph2_0.SetY(std::get<2>(FitParamsPions[i][1]) * 0 + std::get<3>(FitParamsPions[i][1]));
                V_ph2_0.SetZ(0);

                V_ph2_1.SetX(std::get<0>(FitParamsPions[i][1]) * 1 + std::get<1>(FitParamsPions[i][1]));
                V_ph2_1.SetY(std::get<2>(FitParamsPions[i][1]) * 1 + std::get<3>(FitParamsPions[i][1]));
                V_ph2_1.SetZ(1);

                TVector3 V_ph2=V_ph2_1 -V_ph2_0;
                TVector3 Ph2_direction = V_ph2.Unit();


                //True Direction reconstructed energy

                EcalTree->GetEntry(i);

                Ph1_dir_true = Cevent->MomentumPh1().Unit();    //get true momentum from Tree
                Ph2_dir_true = Cevent->MomentumPh2().Unit();

                L_ph1_true.SetPxPyPzE(Ph1_dir_true.X()*EnergyPhoton1[i],
                                      Ph1_dir_true.Y()*EnergyPhoton1[i],
                                      Ph1_dir_true.Z()*EnergyPhoton1[i],
                                      EnergyPhoton1[i]);

                L_ph2_true.SetPxPyPzE(Ph2_dir_true.X()*EnergyPhoton2[i],
                                      Ph2_dir_true.Y()*EnergyPhoton2[i],
                                      Ph2_dir_true.Z()*EnergyPhoton2[i],
                                      EnergyPhoton2[i]);


                //calculate Invariant Mass
                L_pi_true=L_ph1_true + L_ph2_true;
                Double_t invMassPion_true =  L_pi_true.M();
                InvMassSim->Fill(invMassPion_true);

                //true energy reconstructed direction


                L_ph1.SetPxPyPzE(Ph1_direction.X() * Cevent->EnergyPhoton1(),
                                 Ph1_direction.Y() * Cevent->EnergyPhoton1(),
                                 Ph1_direction.Z() * Cevent->EnergyPhoton1(),
                                 Cevent->EnergyPhoton1());
                L_ph2.SetPxPyPzE(Ph2_direction.X() * Cevent->EnergyPhoton2(),
                                 Ph2_direction.Y() * Cevent->EnergyPhoton2(),
                                 Ph2_direction.Z() * Cevent->EnergyPhoton2(),
                                 Cevent->EnergyPhoton2());


                L_pi=L_ph1 + L_ph2;
                Double_t invMassPion=L_pi.M();

                InvMassReco->Fill(invMassPion);

                InvariantMass.push_back(invMassPion);

                //reconstruchted energy reconstruchted direction

                L_ph1_reco.SetPxPyPzE(Ph1_direction.X() * EnergyPhoton1[i],
                                      Ph1_direction.Y() * EnergyPhoton1[i],
                                      Ph1_direction.Z() * EnergyPhoton1[i],
                                      EnergyPhoton1[i]);
                L_ph2_reco.SetPxPyPzE(Ph2_direction.X() * EnergyPhoton2[i],
                                      Ph2_direction.Y() * EnergyPhoton2[i],
                                      Ph2_direction.Z() * EnergyPhoton2[i],
                                      EnergyPhoton2[i]);

                L_pi_reco=L_ph1_reco+L_ph2_reco;

                Double_t invMassPion_reco=L_pi_reco.M();

                InvMassAllreco->Fill(invMassPion_reco);

                //reconstruchted direction

                Double_t anglereco=(360/(2*TMath::Pi()))*Ph1_direction.Angle(Ph2_direction);
                Angle_vs_massReco->Fill(invMassPion, anglereco);

                //simulated direction

                Double_t anglesim=(360/(2*TMath::Pi()))*Ph1_dir_true.Angle(Ph2_dir_true);
                Angle_vs_massSim->Fill(invMassPion_true, anglesim);

                //all reconstructed

                Angle_vs_massAllreco->Fill(invMassPion_reco, anglereco);

                // Fill Histograms for difference between reconstructed and true direction
                delX1->Fill(Ph1_dir_true.X()-Ph1_direction.X());
                delY1->Fill(Ph1_dir_true.Y()-Ph1_direction.Y());
                delZ1->Fill(Ph1_dir_true.Z()-Ph1_direction.Z());

                delX2->Fill(Ph2_dir_true.X()-Ph2_direction.X());
                delY2->Fill(Ph2_dir_true.Y()-Ph2_direction.Y());
                delZ2->Fill(Ph2_dir_true.Z()-Ph2_direction.Z());

        }
        IM1->cd(1);
        InvMassReco->Draw();

        IM1->cd(2);
        InvMassSim->Draw();

        IM1->cd(3);
        InvMassAllreco->Draw();

        IM1->cd(4);
        Angle_vs_massReco->Draw("colz");

        IM1->cd(5);
        Angle_vs_massSim->Draw("colz");

        IM1->cd(6);
        Angle_vs_massAllreco->Draw("colz");


        DirDelta1->cd(1);
        delX1->Draw();

        DirDelta1->cd(3);
        delY1->Draw();

        DirDelta1->cd(5);
        delZ1->Draw();

        DirDelta1->cd(2);
        delX2->Draw();

        DirDelta1->cd(4);
        delY2->Draw();

        DirDelta1->cd(6);
        delZ2->Draw();



}

//------------------------------------------------------------------------------------------------------------------------------------------

// void TROOTAnalysis::MinimizeClosestApproach(){
//
//
//         TCanvas * PiLoc = new TCanvas("MinimizeClosestApproachC");
//         PiLoc->Divide(3,2,0.01,0.01);
//
//         TH1D * pilocX = new TH1D("MinimizeClosestApproachX", "MinimizeClosestApproachX",2000,-500,500);
//         TH1D * pilocY = new TH1D("MinimizeClosestApproachY", "MinimizeClosestApproachY",2000,-500,500);
//         TH1D * pilocZ = new TH1D("MinimizeClosestApproachZ", "MinimizeClosestApproachZ",2000,-500,500);
//
//         TH1D * truelocX = new TH1D("TrueLocationX", "TrueLocationX",2000,-500,500);
//         TH1D * truelocY = new TH1D("TrueLocationY", "TrueLocationY",2000,-500,500);
//         TH1D * truelocZ = new TH1D("TrueLocationZ", "TrueLocationZ",2000,-500,500);
//
//
//
//
//
//
//         Fcn myfcnLoc;
//
//         myfcnLoc.SetParamsPions(FitParamsPions);
//         myfcnLoc.SetMode("approach");
//
//
//
//         TVector3 V_ph1_0_loc;
//         TVector3 V_ph1_1_loc;
//         TVector3 V_ph2_0_loc;
//         TVector3 V_ph2_1_loc;
//
//         for(Int_t i=0; i<nofEntries; i++) {
//
//                 MnUserParameters uparLoc;
//
//                 double error_minimizer_parameters = 1e-4;
//
//                 uparLoc.Add("x", std::get<0>(ClosestApproach[i]), error_minimizer_parameters);
//                 uparLoc.Add("y", std::get<1>(ClosestApproach[i]), error_minimizer_parameters);
//                 uparLoc.Add("z", std::get<2>(ClosestApproach[i]), error_minimizer_parameters);
//                 //
//                 // uparLoc.Add("x", Cevent->GunPos().X(), error_minimizer_parameters);
//                 // uparLoc.Add("y", Cevent->GunPos().X(), error_minimizer_parameters);
//                 //
//                 // uparLoc.Add("z", Cevent->GunPos().X(), error_minimizer_parameters);
//
//
//
//                 EcalTree->GetEntry(i);
//                 std::cout<<"Event: "<< i <<std::endl;
//                 TVector3 Ph1_dir_true_loc=Cevent->MomentumPh1().Unit();
//                 TVector3 Ph2_dir_true_loc=Cevent->MomentumPh2().Unit();
//
//
//
//                 V_ph1_0_loc.SetX(std::get<0>(FitParamsPions[i][0]) * 0 + std::get<1>(FitParamsPions[i][0]));
//                 V_ph1_0_loc.SetY(std::get<2>(FitParamsPions[i][0]) * 0 + std::get<3>(FitParamsPions[i][0]));
//                 V_ph1_0_loc.SetZ(0);
//
//                 V_ph1_1_loc.SetX(std::get<0>(FitParamsPions[i][0]) * 1 + std::get<1>(FitParamsPions[i][0]));
//                 V_ph1_1_loc.SetY(std::get<2>(FitParamsPions[i][0]) * 1 + std::get<3>(FitParamsPions[i][0]));
//                 V_ph1_1_loc.SetZ(1);
//
//                 TVector3 V_ph1_loc=V_ph1_1_loc - V_ph1_0_loc;
//
//                 TVector3 Ph1_direction_loc = V_ph1_loc.Unit();
//
//                 TLorentzVector L_ph1_loc;
//
//                 L_ph1_loc.SetPxPyPzE(Ph1_direction_loc.X() * EnergyPhoton1[i],
//                                      Ph1_direction_loc.Y() * EnergyPhoton1[i],
//                                      Ph1_direction_loc.Z() * EnergyPhoton1[i],
//                                      EnergyPhoton1[i]);
//
//
//                 V_ph2_0_loc.SetX(std::get<0>(FitParamsPions[i][1]) * 0 + std::get<1>(FitParamsPions[i][1]));
//                 V_ph2_0_loc.SetY(std::get<2>(FitParamsPions[i][1]) * 0 + std::get<3>(FitParamsPions[i][1]));
//                 V_ph2_0_loc.SetZ(0);
//
//                 V_ph2_1_loc.SetX(std::get<0>(FitParamsPions[i][1]) * 1 + std::get<1>(FitParamsPions[i][1]));
//                 V_ph2_1_loc.SetY(std::get<2>(FitParamsPions[i][1]) * 1 + std::get<3>(FitParamsPions[i][1]));
//                 V_ph2_1_loc.SetZ(1);
//
//                 TVector3 V_ph2_loc=V_ph2_1_loc -V_ph2_0_loc;
//
//                 TVector3 Ph2_direction_loc = V_ph2_loc.Unit();
//
//                 TLorentzVector L_ph2_loc;
//
//                 L_ph2_loc.SetPxPyPzE(Ph2_direction_loc.X() * EnergyPhoton2[i],
//                                      Ph2_direction_loc.Y() * EnergyPhoton2[i],
//                                      Ph2_direction_loc.Z() * EnergyPhoton2[i],
//                                      EnergyPhoton2[i]);
//
//                 TVector3 startvalues1(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton1[i]),
//                                       std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton1[i]),
//                                       std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton1[i]));
//                 Double_t deltaT1=startvalues1.Mag();
//
//                 TVector3 startvalues2(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton2[i]),
//                                       std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton2[i]),
//                                       std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton2[i]));
//                 Double_t deltaT2=startvalues2.Mag();
//
//
//                 TVector3 PointPh1(std::get<0>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0]),
//                                   std::get<2>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]),
//                                   std::get<2>(ClosestApproach[i]));
//
//                 TVector3 PointPh2(std::get<0>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1]),
//                                   std::get<2>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1]),
//                                   std::get<2>(ClosestApproach[i]));
//
//
//                 myfcnLoc.SetPhotonPoints(PointPh1.X(), PointPh1.Y(), PointPh1.Z(), EnergyPhoton1[i],
//                                          PointPh2.X(), PointPh2.Y(), PointPh2.Z(), EnergyPhoton2[i]);
//                 myfcnLoc.SetDeltas(deltaT1, deltaT2);
//
//                 myfcnLoc.SetshowerCOG(showerCOGPhoton1[i], showerCOGPhoton2[i]);
//
//                 MnMigrad migradLoc(myfcnLoc, uparLoc, 2);
//
//                 FunctionMinimum minLoc = migradLoc();
//
//                 std::cout<<minLoc<<std::endl;
//
//                 MnUserParameterState userParameterStateLoc = minLoc.UserState();
//
//                 pilocX->Fill(userParameterStateLoc.Value("x"));
//                 pilocY->Fill(userParameterStateLoc.Value("y"));
//                 pilocZ->Fill(userParameterStateLoc.Value("z"));
//
//                 truelocX->Fill(std::get<0>(ClosestApproach[i]));
//                 truelocY->Fill(std::get<1>(ClosestApproach[i]));
//                 truelocZ->Fill(std::get<2>(ClosestApproach[i]));
//
//
//
//
//
//
//
//
//
//
//
//         }
//
//
//         PiLoc->cd(1);
//         pilocX->Draw();
//
//         PiLoc->cd(2);
//         pilocY->Draw();
//
//         PiLoc->cd(3);
//         pilocZ->Draw();
//
//         PiLoc->cd(4);
//         truelocX->Draw();
//
//         PiLoc->cd(5);
//         truelocY->Draw();
//
//         PiLoc->cd(6);
//         truelocZ->Draw();
//
// }


void TROOTAnalysis::PionLocator(){



        TCanvas * PiLoc = new TCanvas("MinimizeClosestApproachC");
        PiLoc->Divide(3,3,0.01,0.01);

        TH1D * pilocX = new TH1D("MinimizeClosestApproachX", "MinimizeClosestApproachX",2000,-500,500);
        TH1D * pilocY = new TH1D("MinimizeClosestApproachY", "MinimizeClosestApproachY",2000,-500,500);
        TH1D * pilocZ = new TH1D("MinimizeClosestApproachZ", "MinimizeClosestApproachZ",3000,-2000,500);

        TH1D * inputX = new TH1D("InputLocationX", "InputLocationX",2000,-500,500);
        TH1D * inputY = new TH1D("InputLocationY", "InputLocationY",2000,-500,500);
        TH1D * inputZ = new TH1D("InputLocationZ", "InputLocationZ",3000,-2000,500);

        TH1D * invmasFitX= new TH1D("FitInvariantMassX","FitInvariantMassX", 2000,-500,500);
        TH1D * invmasFitY= new TH1D("FitInvariantMassY","FitInvariantMassY", 2000,-500,500);
        TH1D * invmasFitZ= new TH1D("FitInvariantMassZ","FitInvariantMassZ", 3000,-2000,500);







        Fcn myfcnLoc;
        Fcn myfcnEn;
        myfcnLoc.SetParamsPions(FitParamsPions);

        myfcnLoc.SetMode("approach");

        myfcnEn.SetParamsPions(FitParamsPions);
        myfcnEn.SetMode("InvMass2");

        TVector3 V_ph1_0_loc;
        TVector3 V_ph1_1_loc;
        TVector3 V_ph2_0_loc;
        TVector3 V_ph2_1_loc;

        for(Int_t i=0; i<nofEntries; i++) {

                MnUserParameters uparLoc;
                MnUserParameters uparEn;

                double error_minimizer_parameters = 1e-4;

                uparLoc.Add("x", std::get<0>(ClosestApproach[i]), error_minimizer_parameters);
                uparLoc.Add("y", std::get<1>(ClosestApproach[i]), error_minimizer_parameters);
                uparLoc.Add("z", std::get<2>(ClosestApproach[i]), error_minimizer_parameters);


                //
                // uparLoc.Add("x", Cevent->GunPos().X(), error_minimizer_parameters);
                // uparLoc.Add("y", Cevent->GunPos().X(), error_minimizer_parameters);
                //
                // uparLoc.Add("z", Cevent->GunPos().X(), error_minimizer_parameters);



                EcalTree->GetEntry(i);
                std::cout<<"Event: "<< i <<std::endl;
                TVector3 Ph1_dir_true_loc=Cevent->MomentumPh1().Unit();
                TVector3 Ph2_dir_true_loc=Cevent->MomentumPh2().Unit();



                V_ph1_0_loc.SetX(std::get<0>(FitParamsPions[i][0]) * 0 + std::get<1>(FitParamsPions[i][0]));
                V_ph1_0_loc.SetY(std::get<2>(FitParamsPions[i][0]) * 0 + std::get<3>(FitParamsPions[i][0]));
                V_ph1_0_loc.SetZ(0);

                V_ph1_1_loc.SetX(std::get<0>(FitParamsPions[i][0]) * 1 + std::get<1>(FitParamsPions[i][0]));
                V_ph1_1_loc.SetY(std::get<2>(FitParamsPions[i][0]) * 1 + std::get<3>(FitParamsPions[i][0]));
                V_ph1_1_loc.SetZ(1);

                TVector3 V_ph1_loc=V_ph1_1_loc - V_ph1_0_loc;

                TVector3 Ph1_direction_loc = V_ph1_loc.Unit();

                TLorentzVector L_ph1_loc;

                L_ph1_loc.SetPxPyPzE(Ph1_direction_loc.X() * EnergyPhoton1[i],
                                     Ph1_direction_loc.Y() * EnergyPhoton1[i],
                                     Ph1_direction_loc.Z() * EnergyPhoton1[i],
                                     EnergyPhoton1[i]);


                V_ph2_0_loc.SetX(std::get<0>(FitParamsPions[i][1]) * 0 + std::get<1>(FitParamsPions[i][1]));
                V_ph2_0_loc.SetY(std::get<2>(FitParamsPions[i][1]) * 0 + std::get<3>(FitParamsPions[i][1]));
                V_ph2_0_loc.SetZ(0);

                V_ph2_1_loc.SetX(std::get<0>(FitParamsPions[i][1]) * 1 + std::get<1>(FitParamsPions[i][1]));
                V_ph2_1_loc.SetY(std::get<2>(FitParamsPions[i][1]) * 1 + std::get<3>(FitParamsPions[i][1]));
                V_ph2_1_loc.SetZ(1);

                TVector3 V_ph2_loc=V_ph2_1_loc -V_ph2_0_loc;

                TVector3 Ph2_direction_loc = V_ph2_loc.Unit();

                TLorentzVector L_ph2_loc;

                L_ph2_loc.SetPxPyPzE(Ph2_direction_loc.X() * EnergyPhoton2[i],
                                     Ph2_direction_loc.Y() * EnergyPhoton2[i],
                                     Ph2_direction_loc.Z() * EnergyPhoton2[i],
                                     EnergyPhoton2[i]);

                TVector3 startvalues1(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton1[i]),
                                      std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton1[i]),
                                      std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton1[i]));
                Double_t deltaT1=startvalues1.Mag();

                TVector3 startvalues2(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton2[i]),
                                      std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton2[i]),
                                      std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton2[i]));
                Double_t deltaT2=startvalues2.Mag();


                TVector3 PointPh1(std::get<0>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0]),
                                  std::get<2>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]),
                                  std::get<2>(ClosestApproach[i]));

                TVector3 PointPh2(std::get<0>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1]),
                                  std::get<2>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1]),
                                  std::get<2>(ClosestApproach[i]));


                myfcnLoc.SetPhotonPoints(PointPh1.X(), PointPh1.Y(), PointPh1.Z(), EnergyPhoton1[i],
                                         PointPh2.X(), PointPh2.Y(), PointPh2.Z(), EnergyPhoton2[i]);
                myfcnLoc.SetDeltas(deltaT1, deltaT2);

                myfcnLoc.SetshowerCOG(showerCOGPhoton1[i], showerCOGPhoton2[i]);

                MnMigrad migradLoc(myfcnLoc, uparLoc, 2);

                FunctionMinimum minLoc = migradLoc();

                //std::cout<<minLoc<<std::endl;

                MnUserParameterState userParameterStateLoc = minLoc.UserState();



                inputX->Fill(std::get<0>(ClosestApproach[i]));
                inputY->Fill(std::get<1>(ClosestApproach[i]));
                inputZ->Fill(std::get<2>(ClosestApproach[i]));

                pilocX->Fill(userParameterStateLoc.Value("x"));
                pilocY->Fill(userParameterStateLoc.Value("y"));
                pilocZ->Fill(userParameterStateLoc.Value("z"));

                Double_t error_minimizer_parameters_En=0.0001;

                uparEn.Add("x", userParameterStateLoc.Value("x"), error_minimizer_parameters_En);
                uparEn.Add("y", userParameterStateLoc.Value("y"), error_minimizer_parameters_En);
                uparEn.Add("z", userParameterStateLoc.Value("z"), error_minimizer_parameters_En);

                uparEn.Add("cx1", std::get<0>(showerCOGPhoton1[i]), error_minimizer_parameters_En);
                uparEn.Add("cy1", std::get<1>(showerCOGPhoton1[i]), error_minimizer_parameters_En);
                uparEn.Add("cz1", std::get<2>(showerCOGPhoton1[i]), error_minimizer_parameters_En);

                uparEn.Add("cx2", std::get<0>(showerCOGPhoton2[i]), error_minimizer_parameters_En);
                uparEn.Add("cy2", std::get<1>(showerCOGPhoton2[i]), error_minimizer_parameters_En);
                uparEn.Add("cz2", std::get<2>(showerCOGPhoton2[i]), error_minimizer_parameters_En);

                //uparEn.SetLimits("x",userParameterStateLoc.Value("x")-20, userParameterStateLoc.Value("x")+20);
                //uparEn.SetLimits("x",userParameterStateLoc.Value("y")-30, userParameterStateLoc.Value("y")+30);

                myfcnEn.SetDeltas(deltaT1, deltaT2);
                myfcnEn.SetPhotonPoints(PointPh1.X(), PointPh1.Y(), PointPh1.Z(), EnergyPhoton1[i],
                                        PointPh2.X(), PointPh2.Y(), PointPh2.Z(), EnergyPhoton2[i]);
                myfcnEn.SetshowerCOG(showerCOGPhoton1[i], showerCOGPhoton2[i]);
                myfcnEn.SetClosestApproach(userParameterStateLoc.Value("x"),userParameterStateLoc.Value("y"),userParameterStateLoc.Value("z"));

                MnMigrad migradEn(myfcnEn, uparEn, 2);

                FunctionMinimum minEn = migradEn();

                MnUserParameterState userParameterStateEn=minEn.UserState();


                invmasFitX->Fill(userParameterStateEn.Value("x"));
                invmasFitY->Fill(userParameterStateEn.Value("y"));
                invmasFitZ->Fill(userParameterStateEn.Value("z"));





        }


        PiLoc->cd(4);
        pilocX->Draw();

        PiLoc->cd(5);
        pilocY->Draw();

        PiLoc->cd(6);
        pilocZ->Draw();

        PiLoc->cd(1);
        inputX->Draw();

        PiLoc->cd(2);
        inputY->Draw();

        PiLoc->cd(3);
        inputZ->Draw();

        PiLoc->cd(7);
        invmasFitX->Draw();

        PiLoc->cd(8);
        invmasFitY->Draw();

        PiLoc->cd(9);
        invmasFitZ->Draw();

        // TCanvas * sgd = new TCanvas("sgdC", "sgdC");
        // TH1D * lkl = new TH1D("lkl", "lkl", 1000,0,1000);
        //
        // Double_t x[1000];
        // Double_t y[1000];
        // for(Int_t j=0; j<nofEntries; j++) {
        //
        //         EcalTree->GetEntry(j);
        //
        //         TVector3 origdir1=Cevent->MomentumPh1().Unit();
        //         TVector3 origdir2=Cevent->MomentumPh2().Unit();
        //
        //
        //
        //         TVector3 dir_ph1(std::get<0>(showerCOGPhoton1[j])-std::get<0>(ClosestApproach[j]),
        //                          std::get<1>(showerCOGPhoton1[j])-std::get<1>(ClosestApproach[j]),
        //                          std::get<2>(showerCOGPhoton1[j])-std::get<2>(ClosestApproach[j]));
        //         dir_ph1=dir_ph1.Unit();
        //
        //         TLorentzVector lv_ph1;
        //         lv_ph1.SetPxPyPzE(dir_ph1.X()*Cevent->EnergyPhoton1(),dir_ph1.Y()*Cevent->EnergyPhoton1(),dir_ph1.Z()*Cevent->EnergyPhoton1(),Cevent->EnergyPhoton1());
        //
        //         TVector3 dir_ph2(std::get<0>(showerCOGPhoton2[j])-std::get<0>(ClosestApproach[j]),
        //                          std::get<1>(showerCOGPhoton2[j])-std::get<1>(ClosestApproach[j]),
        //                          std::get<2>(showerCOGPhoton2[j])-std::get<2>(ClosestApproach[j]));
        //         dir_ph2=dir_ph2.Unit();
        //
        //         TLorentzVector lv_ph2;
        //         lv_ph2.SetPxPyPzE(dir_ph2.X()*Cevent->EnergyPhoton2(),dir_ph2.Y()*Cevent->EnergyPhoton2(),dir_ph2.Z()*Cevent->EnergyPhoton2(),Cevent->EnergyPhoton2());
        //
        //         lv_ph1.Print();
        //         lv_ph2.Print();
        //         TLorentzVector lv_pi;
        //         lv_pi  = lv_ph1 + lv_ph2;
        //
        //         lv_pi.Print();
        //
        //         Double_t curr_invmass=lv_pi.M();
        //         std::cout<<curr_invmass<<std::endl;
        //         lkl->Fill(curr_invmass);
        //
        //
        //         // origdir1.Print();
        //         // dir_ph1.Print();
        //         //
        //         // std::cout<<"PH2____________________________"<<std::endl;
        //         //
        //         // origdir2.Print();
        //         // dir_ph2.Print();
        //
        //         std::cout<<"PH1____________________________"<<std::endl;
        //
        // }
        // // TGraph* g1= new TGraph(1000, x,y);
        // sgd->cd(0);
        // lkl->Draw();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

void TROOTAnalysis::PlotChimap(Int_t event){

        TCanvas * ChiMap=new TCanvas("ChiMapC","Chimap");
        ChiMap->Divide(2,1,0.01,0.01);
        TH3D * chimap = new TH3D("chimap", "chimap",1000,-500,500,1000,-500,500,1000,-1600,400);



        for(Int_t i=0; i<nofEntries; i++) {

                EcalTree->GetEntry(i);

                // Double_t e1=Cevent->EnergyPhoton1();
                // Double_t e2=Cevent->EnergyPhoton2();
                Double_t e1=EnergyPhoton1[i];
                Double_t e2=EnergyPhoton2[i];
                TVector3 startvalues1(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton1[i]),
                                      std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton1[i]),
                                      std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton1[i]));
                Double_t deltaT1=startvalues1.Mag();

                TVector3 startvalues2(std::get<0>(ClosestApproach[i]) - std::get<0>(showerCOGPhoton2[i]),
                                      std::get<1>(ClosestApproach[i]) - std::get<1>(showerCOGPhoton2[i]),
                                      std::get<2>(ClosestApproach[i]) - std::get<2>(showerCOGPhoton2[i]));
                Double_t deltaT2=startvalues2.Mag();

                TVector3 PointPh1(std::get<0>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][0]),
                                  std::get<2>(FitParamsPions[i][0]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][0]),
                                  std::get<2>(ClosestApproach[i]));

                TVector3 PointPh2(std::get<0>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<1>(FitParamsPions[i][1]),
                                  std::get<2>(FitParamsPions[i][1]) * std::get<2>(ClosestApproach[i]) + std::get<3>(FitParamsPions[i][1]),
                                  std::get<2>(ClosestApproach[i]));

                PointPh1.Print();
                PointPh2.Print();
                std::cout<<std::get<0>(ClosestApproach[i])<<":"<<std::get<1>(ClosestApproach[i])<<":"<<std::get<2>(ClosestApproach[i])<<std::endl;

                Double_t chisq;

                for(Int_t x = -500; x<501; x+=2) {
                        for(Int_t y = -500; y<501; y+=2) {
                                for(Int_t z = -1600; z< 401; z+=2) {

                                        Double_t sigPh1=(1/TMath::Sqrt(e1))*deltaT1;
                                        Double_t sigPh2=(1/TMath::Sqrt(e2))*deltaT2;


                                        Double_t delR1=TMath::Sqrt(((PointPh1.X()-x)*(PointPh1.X()-x)) + ((PointPh1.Y()-y)*(PointPh1.Y()-y)) + ((PointPh1.Z()-z)*(PointPh1.Z()-z)));
                                        Double_t delR2=TMath::Sqrt(((PointPh2.X()-x)*(PointPh2.X()-x)) + ((PointPh2.Y()-y)*(PointPh2.Y()-y)) + ((PointPh2.Z()-z)*(PointPh2.Z()-z)));

                                        TVector3 Dir_ph1(x-std::get<0>(showerCOGPhoton1[i]),y-std::get<1>(showerCOGPhoton1[i]),z-std::get<2>(showerCOGPhoton1[i]));

                                        Dir_ph1=Dir_ph1.Unit();

                                        TLorentzVector Lv_ph1(Dir_ph1.X()*e1,Dir_ph1.Y()*e1,Dir_ph1.Z()*e1,e1);


                                        TVector3 Dir_ph2(x-std::get<0>(showerCOGPhoton2[i]),y-std::get<1>(showerCOGPhoton2[i]),z-std::get<2>(showerCOGPhoton2[i]));

                                        Dir_ph2=Dir_ph2.Unit();

                                        TLorentzVector Lv_ph2(Dir_ph2.X()*e2,Dir_ph2.Y()*e2,Dir_ph2.Z()*e2,e2);

                                        TLorentzVector Lv_pi=Lv_ph1+Lv_ph2;

                                        Double_t Curr_invmass=Lv_pi.M();




                                        chisq=((delR1/sigPh1)*(delR1/sigPh1))+((delR2/sigPh2)*(delR2/sigPh2))+(((Curr_invmass-134.9766)/5)*((Curr_invmass-134.9766)/5));

                                        chimap->Fill(x,y,z,chisq);


                                }
                                //std::cout<<y<<std::endl;
                        }

                }
                std::cout<<"projection"<<std::endl;
                ChiMap->cd(1);
                chimap->Project3D("xz")->Draw("colz");
                std::cout<<"projection"<<std::endl;
                ChiMap->cd(2);
                chimap->Project3D("yz")->Draw("colz");

                std::string ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy())+"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".pdf";

                ChiMap->Print(ChiMapPath.c_str());

                ChiMapPath="Chi2Mapz/ChiMap" +std::to_string(Cevent->GapEnergy()) +"_" +std::to_string(Cevent->GunPos().Z())+"_"+ std::to_string(i) + ".C";

                ChiMap->Print(ChiMapPath.c_str());

                chimap->Reset();




        }

}



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



void TROOTAnalysis::PrintFitHists(Int_t minevent, Int_t maxevent){
        Int_t ent=nofEntries;
        TCanvas * c1 = new TCanvas("Fit", "Fit", 1500,1300);

        c1->Divide(2,2,0.01, 0.01);

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
        gStyle->SetOptStat(2222);
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

        //c1->cd(5);
        //er1->Draw();

        //c1->cd(6);
        //er2->Draw();

        if(pathset) {
                //
                std::string imgname= "Fit.png";
                std::string imgpath = savepath + '/' + imgname;
                // TImage * img2 =TImage::Create();
                // img2->FromPad(c1);
                // img2->WriteImage(imgpath.c_str());
                // delete img2;

                c1->Print(imgpath.c_str());

                imgname= "Fit.C";
                imgpath = savepath + '/' + imgname;

                c1->Print(imgpath.c_str());
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

        for(Int_t eventstodo= 0; eventstodo<nofEntries; eventstodo++) {    //loop over all simulated events

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




        //transform from copynumber coordinates to geant4 coordinates
        std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > Transcoglist;
        std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t, Double_t> > > TransfomedCOGs;




        Double_t offsetX=((-1*calsizeXY/2)+tiledimX/2);
        Double_t offsetY=((-1*calsizeXY/2)+tiledimY/2);
        Double_t offsetZ=((((AbsoThickness+GapThickness)*(-1)*nofLayers)/2)+(AbsoThickness+GapThickness/2));
        std::cout<<"X: "<<offsetX<<"Y:"<<offsetY<<"Z: "<<offsetZ<<std::endl;

        for(Int_t i = 0; i<nofEntries; i++) {                                         //transformation to geant4 coordinate system
                for(Int_t j = 0; j<COGCollection[i].size(); j++) {

                        auto tc = std::make_tuple((offsetX + std::get<0>(COGCollection[i][j]) * tiledimX),                // X
                                                  (offsetY + std::get<1>(COGCollection[i][j]) * tiledimY),                 // Y
                                                  (offsetZ + std::get<2>(COGCollection[i][j]) * (AbsoThickness+GapThickness)),                               // Z
                                                  std::get<3>(COGCollection[i][j]),
                                                  std::get<4>(COGCollection[i][j]),
                                                  std::get<5>(COGCollection[i][j])  );

                        Transcoglist.push_back(tc);

                }
                TransfomedCOGs.push_back(Transcoglist);
                Transcoglist.clear();
        }





        Fcn myfcn;
        myfcn.SetCOGs(TransfomedCOGs);
        myfcn.SetMode("photon");
        //myfcn.PrintCOGs();
        MnUserParameters upar;
        double error_minimizer_parameters = 1e-4;

        upar.Add("mx", 0., error_minimizer_parameters);
        upar.Add("tx", 1., error_minimizer_parameters);
        upar.Add("my", 0., error_minimizer_parameters);
        upar.Add("ty", 1., error_minimizer_parameters);

        cout << upar << endl;

        MnMigrad migrad(myfcn, upar, 2);



        for(int events = 0; events < nofEntries; events++) {

                upar.SetValue("mx", 0.);
                upar.SetValue("tx", 0.);
                upar.SetValue("my", 0.);
                upar.SetValue("ty", 0.);

                myfcn.SetCurrentEvent(events);

                std::cout<<"event: "<<events<<std::endl;
                FunctionMinimum min = migrad();

                std::cout<<min<<std::endl;


                //MnUserCovariance cov = min.UserCovariance();

                //std::cout<<cov(0 ,0)<<" "<<cov(1,1)<<" "<<cov(2,2)<<" "<<cov(3,3)<<std::endl;


                MnUserParameterState userParameterState = min.UserState();

                auto tp = std::make_tuple(userParameterState.Value("mx"), userParameterState.Value("tx"),
                                          userParameterState.Value("my"), userParameterState.Value("ty"));
                FitParamsGamma.push_back(tp);



                // Double_t covar[4][4];
                // Double_t error[4];
                //
                // error[0]= userParameterState.Error("mx");
                // error[1]= userParameterState.Error("tx");
                // error[2]= userParameterState.Error("my");
                // error[3]= userParameterState.Error("ty");
                //
                // //Fill correlation Histograms
                // if(COGCollection[events].size() > 1) {
                //
                //
                //         for(Int_t row=0; row<4; row++) {
                //                 for(Int_t col=0; col<4; col++) {
                //                         covar[row][col]=cov(row,col);
                //
                //                 }
                //         }
                //
                //         co1->Fill(covar[0][1]/(error[0]*error[1]));
                //         co2->Fill(covar[0][2]/(error[0]*error[2]));
                //         co3->Fill(covar[0][3]/(error[0]*error[3]));
                //         co4->Fill(covar[1][2]/(error[1]*error[2]));
                //         co5->Fill(covar[1][3]/(error[1]*error[3]));
                //         co6->Fill(covar[2][3]/(error[2]*error[3]));
                // }
                // auto tp=std::make_tuple(userParameterState.Value("a1"),userParameterState.Value("a2"),userParameterState.Value("a3"),
                //                         userParameterState.Value("x1"),userParameterState.Value("x2"),userParameterState.Value("x3"));
                //
                // FitParams.push_back(tp);
        }



        // corr1->cd(1);
        // co1->Draw("");
        //
        // corr1->cd(2);
        // co2->Draw("");
        //
        // corr1->cd(3);
        // co3->Draw("");
        //
        // corr1->cd(4);
        // co4->Draw("");
        //
        // corr1->cd(5);
        // co5->Draw("");
        //
        // corr1->cd(6);
        // co6->Draw("");

}
