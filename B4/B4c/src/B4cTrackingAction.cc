#include "B4cTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "B4cTrackInformation.hh"

B4cTrackingAction::B4cTrackingAction() : G4UserTrackingAction()
{
        ;
}

B4cTrackingAction::~B4cTrackingAction()
{
        ;
}

void B4cTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
        if(aTrack->GetTrackID()==2)
        {
                isdecayed=true;
                G4cout<<"TrackID1: "<<aTrack->GetTrackID()<<G4endl;
                //G4cout<<"true1"<<G4endl;
                G4cout<<aTrack->GetVolume()->GetName()<<G4endl;

                B4cTrackInformation* anInfo = new B4cTrackInformation(aTrack);
                anInfo->SetOriginalPhotonNumber(1);
                G4cout<<anInfo->GetOriginalPhotonNumber()<<G4endl;
                G4cout<<"anInfo:"<<G4endl;
                anInfo->Print();
                G4Track* theTrack = (G4Track*)aTrack;
                theTrack->SetUserInformation(anInfo);
                //G4cout<<"true2"<<G4endl;


        }

        if(aTrack->GetTrackID()==3)
        {
                isdecayed=true;
                G4cout<<"TrackID2: "<<aTrack->GetTrackID()<<G4endl;
                //G4cout<<"true1"<<G4endl;
                G4cout<<aTrack->GetVolume()->GetName()<<G4endl;

                B4cTrackInformation* anInfo = new B4cTrackInformation(aTrack);
                anInfo->SetOriginalPhotonNumber(2);
                G4cout<<anInfo->GetOriginalPhotonNumber()<<G4endl;

                G4cout<<"anInfo:"<<G4endl;
                anInfo->Print();
                G4Track* theTrack = (G4Track*)aTrack;
                theTrack->SetUserInformation(anInfo);
                //G4cout<<"true2"<<G4endl;


        }

}

void B4cTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
        //G4cout<<"true3"<<G4endl;
        if(isdecayed) {
                G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
                if(secondaries)
                {
                        //G4cout<<"true4"<<G4endl;
                        B4cTrackInformation* info = (B4cTrackInformation*)(aTrack->GetUserInformation());
                        size_t nSeco = secondaries->size();
                        //G4cout<<nSeco<<G4endl;
                        if(nSeco>0)
                        {
                                //G4cout<<"true5"<<G4endl;
                                for(size_t i=0; i<nSeco; i++)
                                {
                                        B4cTrackInformation* infoNew = new B4cTrackInformation(info);
                                        //infoNew->Print();
                                        (*secondaries)[i]->SetUserInformation(infoNew);
                                }
                        }
                }
        }
}
