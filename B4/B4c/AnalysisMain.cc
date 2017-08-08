// Main function for analyzing the TTree

#include <iostream>
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"
#include "ROOTAnalysis.hh"
#include <stdlib.h>


int main(int argc, char * argv[]) {

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();


//
//    gStyle->SetCanvasBorderMode(0);
//    gStyle->SetPadBorderMode(0);
//
//    gStyle->SetPalette(1);
//    gStyle->SetOptStat(1);
//    gStyle->SetOptFit(11);
//    //gStyle->SetOptTitle(0);
//    gStyle->SetStatBorderSize(1);
//    gStyle->SetStatColor(10);
//    gStyle->SetCanvasColor(10);
//    gStyle->SetPadLeftMargin(0.16);
//    gStyle->SetPadBottomMargin(0.16);
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
//    gStyle->SetOptTitle(0);
//    gStyle->SetTitleSize(0.048,"xy");
//    gStyle->SetLabelSize(0.04,"xy")
//    gStyle->SetTitleOffset(1.3,"x");
//    gStyle->SetTitleOffset(1.3,"y");


        TChain * ch1 = new TChain("eventTree");


        Int_t minl = std::stoi(argv[1]);
        Int_t maxl = std::stoi(argv[2]);
        Int_t minevt = std::stoi(argv[3]);
        Int_t maxevt = std::stoi(argv[4]);


        std::cout << "Adding tree to the Chain"<< std::endl;

        if(argc>6) {
                char * treepath = argv[6];
                ch1->Add(treepath);
                ch1->Draw("");
                TROOTAnalysis A(ch1);
                std::cout<<"Added Tree"<<std::endl;

                std::cout << "Created Analysis Class" << std::endl;

                std::string s(argv[7]);

                A.SetPath(s);

                //A.plotEvent(100);
                //A.findShowercenter(104, 105);
                A.CalcCOG(minl-1, maxl, minevt, maxevt);
                A.FitCOGs(minevt, maxevt, 10 );

                A.PlotProjection(100);
                //A.PrintFitParams();
                A.PrintFitHists(minevt, maxevt);

                //A.CleanCOGs(minl, maxl, minevt, maxevt);
        }
        else if(argc = 6) {
                TApplication* app = new TApplication("app", 0, 0, 0);

                std::cout<<"aereargerhrehg"<<std::endl;

                ch1->Add("Pion.root");

                std::cout<<"Added Tree"<<std::endl;

                ch1->Draw("");

                std::cout<<"oiueroiueio"<<std::endl;

                TROOTAnalysis A(ch1);


                std::cout << "Created Analysis Class" << std::endl;

                //A.plotEvent(100);
                //A.findShowercenter(104, 105);
                //  A.CalcCOG(minl-1, maxl, minevt, maxevt);
                //  A.FitCOGs(minevt, maxevt, 10 );

                //  A.PlotProjection(100);
                //A.PrintFitParams();
                //  A.PrintFitHists(minevt, maxevt);

                //A.CleanCOGs(minl, maxl, minevt, maxevt);
                app->Run();
        }
        else{std::cout<<"error"<<std::endl; }




//A.PrintERes();
//
// if(std::stoi(argv[5])==0){
//
//
// }
//
// else{
//   A.CalcCOGwithFit(minevt, maxevt);
//}

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Computing took "
                  << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
                  <<" seconds"<<std::endl;




        return 0;
}
