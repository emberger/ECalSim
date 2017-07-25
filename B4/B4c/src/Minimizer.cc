#include "Minimizer.hh"



double Fcn::operator()(const std::vector<double> &par) const
{

  double chisq = 0.;
  std::cout<<COGs[currentEvent].size()<<std::endl;

  for(int i = 0; i<COGs[currentEvent].size(); i++){

        double termx = 0;
        double termy = 0;

          termx += ((std::get<0>(COGs[currentEvent][i]) - (par[0] * std::get<2>(COGs[currentEvent][i])+par[1])) / std::get<3>(COGs[currentEvent][i]))
                    *std::get<5>(COGs[currentEvent][i]);

          termy += ((std::get<1>(COGs[currentEvent][i]) - (par[2] * std::get<2>(COGs[currentEvent][i])+par[3])) / std::get<4>(COGs[currentEvent][i]))
                    *std::get<5>(COGs[currentEvent][i]);


          chisq +=termx*termx+termy*termy ;
    }

      // XYZVector xp;
      // xp.SetCoordinates(std::get<0>(COGs[currentEvent][i]),std::get<1>(COGs[currentEvent][1]),std::get<2>(COGs[currentEvent][i]));
      // XYZVector x0(par[0], par[1], par[2]);
      // XYZVector x1(par[3], par[4], par[5]);
      //
      // XYZVector u = x1 .Unit();
      //
      //
      // dist += TMath::Sqrt(((x0-xp)-((x0-xp).Dot(u)) * u ) .Mag2());
  // }


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
