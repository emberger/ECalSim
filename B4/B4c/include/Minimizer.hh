#ifndef MINIMIZER_HH
#define MINIMIZER_HH

#include <iostream>
#include <vector>
// root include
#include <TFile.h>
#include <TF1.h>

//minuit includes
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FCNBase.h"


using namespace ROOT::Minuit2;
using namespace std;



class Fcn : public ROOT::Minuit2::FCNBase
{
public:
  //
  // /** Funktion die minimiert werden soll */
   virtual double operator()(const std::vector<double> &par) const;
  //
  // /** Error definition: Der Fehler e auf einen Parameter x ist definiert
  // * so dass f(x+-e) = f(x) + Up().
  // * Für einen Chi^2 fit ist 1 das richtige um 1 sigma Fehler zu bekommen. */
   virtual double Up() const;
  //
   inline void SetCurrentEvent(Int_t ce){currentEvent=ce;};
   inline void SetCOGs(std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t>>> c){COGs=c;};
   void PrintCOGs();

private:
  std::vector<std::vector<std::tuple<Double_t,Double_t, Double_t, Double_t, Double_t>>> COGs;
  Int_t currentEvent;

};

// int Fit2Sigma(TH1D* hist){
//   hist->Fit("gaus", "Q");
//   double sigma = hist->GetFunction("gaus")->GetParameter(2);
//   hist->Fit("gaus", "Q", "", -2*sigma, 2*sigma);
//   return 0;
// }




#endif  // MINIMIZER_H
