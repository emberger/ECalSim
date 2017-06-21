// Main function for analyzing the TTree

#include <iostream>
#include "TApplication.h"
#include "TStyle.h"
#include "TChain.h"
#include "ROOTAnalysis.hh"
#include <stdlib.h>


int main(int argc, char const *argv[]) {

 TApplication* app = new TApplication("app", 0, 0, 0);
//
//   	gStyle->SetCanvasBorderMode(0);
//   	gStyle->SetPadBorderMode(0);
//
//   	gStyle->SetPalette(1);
//   	gStyle->SetOptStat(1);
//   	gStyle->SetOptFit(11);
//   	//gStyle->SetOptTitle(0);
//   	gStyle->SetStatBorderSize(1);
//   	gStyle->SetStatColor(10);
//   	gStyle->SetCanvasColor(10);
//   	gStyle->SetPadLeftMargin(0.16);
//   	gStyle->SetPadBottomMargin(0.16);
//   	gStyle->SetPadTickX(1);
//   	gStyle->SetPadTickY(1);
//   	gStyle->SetOptTitle(0);
//   	gStyle->SetTitleSize(0.048,"xy");
//   	gStyle->SetLabelSize(0.04,"xy");
//   	gStyle->SetTitleOffset(1.3,"x");
//   	gStyle->SetTitleOffset(1.3,"y");

TChain * ch1 = new TChain("eventTree");

std::cout << "Adding " <<"ECalEventTree.root"<<" to the Chain"<< std::endl;
ch1->Add("ECalEventTree.root");
ch1->Draw("");
std::cout<<"Added Tree"<<std::endl;

TROOTAnalysis A(ch1);
std::cout << "Created Analysis Class" << std::endl;

Int_t evt= std::stoi(argv[1]);
Int_t ml= std::stoi(argv[2]);

//A.plotEvent(evt-1);
A.GetCOG(ml, evt-1);
app->Run();

  return 0;
}
