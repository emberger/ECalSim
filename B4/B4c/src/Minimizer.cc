#include "Minimizer.hh"


double Fcn::operator()(const std::vector<double> &par) const
{

        double chisq = 0.;
        //std::cout<<COGs[currentEvent].size()<<std::endl;
        if(mode=="photon") {
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

        else if(mode=="pion") {

                Double_t X_d=((std::get<0>(ParamsPions[currentEvent][0])*par[0] + std::get<1>(ParamsPions[currentEvent][0]))-
                              (std::get<0>(ParamsPions[currentEvent][1])*par[0] + std::get<1>(ParamsPions[currentEvent][1])));

                Double_t Y_d=((std::get<2>(ParamsPions[currentEvent][0])*par[0] + std::get<3>(ParamsPions[currentEvent][0]))-
                              (std::get<2>(ParamsPions[currentEvent][1])*par[0] + std::get<3>(ParamsPions[currentEvent][1])));

                chisq=TMath::Sqrt( X_d*X_d + Y_d*Y_d);



        }

        else if(mode=="approach") {
                Double_t sigPh1=(1/TMath::Sqrt(l_ph1[3]))*delta1;
                Double_t sigPh2=(1/TMath::Sqrt(l_ph2[3]))*delta2;

                //std::cout<<"sig1: "<<sigPh1<<" sig2: "<<sigPh2<<std::endl;


                Double_t delR1=TMath::Sqrt(((l_ph1[0]-par[0])*(l_ph1[0]-par[0])) + ((l_ph1[1]-par[1])*(l_ph1[1]-par[1])) + ((l_ph1[2]-par[2])*(l_ph1[2]-par[2])));
                Double_t delR2=TMath::Sqrt(((l_ph2[0]-par[0])*(l_ph2[0]-par[0])) + ((l_ph2[1]-par[1])*(l_ph2[1]-par[1])) + ((l_ph2[2]-par[2])*(l_ph2[2]-par[2])));

                chisq=(delR1/sigPh1)*(delR1/sigPh1)+(delR2/sigPh2)*(delR2/sigPh2);

        }

        else if(mode=="InvMass") {
                Double_t sigPh1=(1/TMath::Sqrt(l_ph1[3]))*delta1;
                Double_t sigPh2=(1/TMath::Sqrt(l_ph2[3]))*delta2;

                //std::cout<<"sig1: "<<sigPh1<<" sig2: "<<sigPh2<<std::endl;


                Double_t delR1=TMath::Sqrt(((l_ph1[0]-par[0])*(l_ph1[0]-par[0])) + ((l_ph1[1]-par[1])*(l_ph1[1]-par[1])) + ((l_ph1[2]-par[2])*(l_ph1[2]-par[2])));
                Double_t delR2=TMath::Sqrt(((l_ph2[0]-par[0])*(l_ph2[0]-par[0])) + ((l_ph2[1]-par[1])*(l_ph2[1]-par[1])) + ((l_ph2[2]-par[2])*(l_ph2[2]-par[2])));


                TVector3 dir_ph1(std::get<0>(cogph1)-par[0],std::get<1>(cogph1)-par[1],std::get<2>(cogph1)-par[2]);

                dir_ph1=dir_ph1.Unit();

                TLorentzVector lv_ph1(dir_ph1.X()*l_ph1[3],dir_ph1.Y()*l_ph1[3],dir_ph1.Z()*l_ph1[3],l_ph1[3]);



                TVector3 dir_ph2(std::get<0>(cogph2)-par[0],std::get<1>(cogph2)-par[1],std::get<2>(cogph2)-par[2]);

                dir_ph2=dir_ph2.Unit();

                TLorentzVector lv_ph2(dir_ph2.X()*l_ph2[3],dir_ph2.Y()*l_ph2[3],dir_ph2.Z()*l_ph2[3],l_ph2[3]);

                TLorentzVector lv_pi;
                lv_pi=lv_ph1+lv_ph2;

                Double_t curr_invmass=lv_pi.M();

                std::cout<<par[0]<<":"<<par[1]<<":"<<par[2]<<std::endl;


                chisq=((delR1/sigPh1)*(delR1/sigPh1))+((delR2/sigPh2)*(delR2/sigPh2))+(((curr_invmass-134.9766)/0.001)*((curr_invmass-134.9766)/0.001));


        }

        else if(mode=="InvMass2") {

                TVector3 dir_ph1(par[0]-par[3], par[1]-par[4], par[2]-par[5]);
                dir_ph1=dir_ph1.Unit();

                TLorentzVector lv_ph1(dir_ph1.X()*l_ph1[3],dir_ph1.Y()*l_ph1[3],dir_ph1.Z()*l_ph1[3],l_ph1[3]);

                TVector3 dir_ph2(par[0]-par[6], par[1]-par[7], par[2]-par[8]);
                dir_ph2=dir_ph2.Unit();

                TLorentzVector lv_ph2(dir_ph2.X()*l_ph2[3],dir_ph2.Y()*l_ph2[3],dir_ph2.Z()*l_ph2[3],l_ph2[3]);

                TLorentzVector lv_pi=lv_ph1+lv_ph1;

                Double_t curr_invmass=lv_pi.M();

                std::cout<<par[0]<<":"<<par[1]<<":"<<par[2]<<std::endl;

                chisq=(((curr_invmass-134.9766)/0.1)*((curr_invmass-134.9766)/0.1));


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
