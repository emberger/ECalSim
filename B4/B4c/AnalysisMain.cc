// Main function for analyzing the TTree

#include <iostream>
#include "TApplication.h"
#include "TChain.h"
#include "ROOTAnalysis.hh"



int main(int argc, char const *argv[]) {

TApplication * app = new TApplication("app", 0,0,0);

TChain * ch1 = new TChain("event Tree");
std::cout << "Adding" <<argv[1]<<"to the Chain"<< std::endl;
ch1->Add(argv[1]);



  return 0;
}
