#include "Minimizer.hh"



double Fcn::operator()(const std::vector<double> &par) const
{

        double chisq = 0.;
        //std::cout<<COGs[currentEvent].size()<<std::endl;
        if(approachmode==false) {
                for(int i = 0; i<COGs[currentEvent].size(); i++) {

                        double termx = 0;
                        double termy = 0;

                        termx += ((std::get<0>(COGs[currentEvent][i]) - (par[0] * std::get<2>(COGs[currentEvent][i])+par[1])) / std::get<3>(COGs[currentEvent][i]))
                                 *std::get<5>(COGs[currentEvent][i]);

                        termy += ((std::get<1>(COGs[currentEvent][i]) - (par[2] * std::get<2>(COGs[currentEvent][i])+par[3])) / std::get<4>(COGs[currentEvent][i]))
                                 *std::get<5>(COGs[currentEvent][i]);


                        chisq +=termx*termx+termy*termy;
                }
        }

        else if(approachmode==true) {

                Double_t X_d=((std::get<0>(ParamsPions[currentEvent][0])*par[0] + std::get<1>(ParamsPions[currentEvent][0]))-
                              (std::get<0>(ParamsPions[currentEvent][1])*par[0] + std::get<1>(ParamsPions[currentEvent][1])));

                Double_t Y_d=((std::get<2>(ParamsPions[currentEvent][0])*par[0] + std::get<3>(ParamsPions[currentEvent][0]))-
                              (std::get<2>(ParamsPions[currentEvent][1])*par[0] + std::get<3>(ParamsPions[currentEvent][1])));

                chisq=TMath::Sqrt( X_d*X_d + Y_d*Y_d);



        }



        return chisq;


}

double Fcn::Up() const
{
        return 1;
}


void Fcn::PrintCOGs(){

        for(Int_t p=0; p<10; p++) {
                for(Int_t q=0; q<COGs[p].size(); q++) {
                        std::cout<<std::get<0>(COGs[p][q])<<" : "<<std::get<1>(COGs[p][q])<<" : "<<std::get<2>(COGs[p][q])
                                 <<" ERRORX: "<<std::get<3>(COGs[p][q])<<" ERRORY: "<<std::get<4>(COGs[p][q])<<std::endl;

                }
                std::cout<<"|||||||||||||||||||"<<std::endl;
        }
}
