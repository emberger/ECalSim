#include "Minimizer.hh"



double Fcn::operator()(const std::vector<double> &par) const
{

  double chisq = 0.;

    for(int i = 0; i<COGs[currentEvent].size(); i++){

        double term = 0;

          term += (std::get<0>(COGs[currentEvent][i]) - (par[0] * std::get<2>(COGs[currentEvent][i])+par[1])) / std::get<3>(COGs[currentEvent][i]);

          //term += (std::get<1>(COGs[currentEvent][i]) - (par[2] * std::get<2>(COGs[currentEvent][i])+par[3])) / std::get<4>(COGs[currentEvent][i]);


          chisq += term * term;
        }
        for(int i = 0; i<COGs[currentEvent].size(); i++){

            double term = 0;

              //term += (std::get<0>(COGs[currentEvent][i]) - (par[0] * std::get<2>(COGs[currentEvent][i])+par[1])) / std::get<3>(COGs[currentEvent][i]);

              term += (std::get<1>(COGs[currentEvent][i]) - (par[2] * std::get<2>(COGs[currentEvent][i])+par[3])) / std::get<4>(COGs[currentEvent][i]);


              chisq += term * term;
            }


        return chisq;
      

}

double Fcn::Up() const
{
  return 1;
}


void Fcn::PrintCOGs(){

  for(Int_t p=0;p<10;p++){
    for(Int_t q=0;q<COGs[p].size();q++){
  std::cout<<std::get<0>(COGs[p][q])<<" : "<<std::get<1>(COGs[p][q])<<" : "<<std::get<2>(COGs[p][q])
  <<" ERRORX: "<<std::get<3>(COGs[p][q])<<" ERRORY: "<<std::get<4>(COGs[p][q])<<std::endl;

  }
  std::cout<<"|||||||||||||||||||"<<std::endl;
  }
}
