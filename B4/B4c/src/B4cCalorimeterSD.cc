//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4cCalorimeterSD.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file B4cCalorimeterSD.cc
/// \brief Implementation of the B4cCalorimeterSD class

#include "B4cCalorimeterSD.hh"
#include "B4cReadoutGeometry.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "B4cDetParams.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorimeterSD::B4cCalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorimeterSD::~B4cCalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection 
    = new B4cCalorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

if(this->GetName()=="GapSD"){
	//create Hit
	//Calculate number of collection cells, add additional ones for each layer and another one for total accounting

  G4int nofEntries=GetInst().GetfNofLayers() * (GetInst().GetnofTilesX() * GetInst().GetnofTilesY()) + GetInst().GetfNofLayers() + 1;


  for (G4int i=0; i<nofEntries; i++ ) {
    fHitsCollection->insert(new B4cCalorHit());
  }
}

else{
	// Create hits
	  // fNofCells for cells + one more for total sums

	for (G4int i=0; i<fNofCells+1; i++ ) {
	    fHitsCollection->insert(new B4cCalorHit());
	}
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B4cCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory* ROhist)
{  


//	if(this->GetName()=="GapSD"){
//		auto tmp=this->GetROgeometry();
//		auto tmp2= tmp ->CheckROVolume(step, ROhist);
//
//		G4cout<<this->GetROgeometry()-> GetROWorld()<<"\t"<<this->GetROgeometry()-> GetName()<<"\t"<< this->isActive() <<G4endl;
//		G4cout <<tmp2<<G4endl;
//		G4cout<<ROhist<<G4endl;
//	}


  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( edep==0. && stepLength == 0. ) return false;      

  auto touchable = (step->GetPreStepPoint()->GetTouchable());
    
  // Get calorimeter cell id 
  auto layerNumber = touchable->GetReplicaNumber(1);
  
//  G4cout<<layerNumber<<G4endl;


  if(this->GetName()=="GapSD"){
//
//  //Get copynumbers to specify cell
//
  auto Cell=ROhist->GetReplicaNumber();
  auto CellV=ROhist-> GetVolume()->GetName();

  auto Strip=ROhist->GetReplicaNumber(1);
  auto StripV=ROhist-> GetVolume(1)->GetName();

  auto Gap=ROhist->GetReplicaNumber(2);
  auto GapV=ROhist-> GetVolume(2)->GetName();

  auto Layer=ROhist->GetReplicaNumber(3);
  auto LayerV=ROhist-> GetVolume(3)->GetName();


  //Calculate CellID

  G4int ROCellID=Layer*GetInt.GettilesPerLayer()


//  auto stageFour=ROhist->GetReplicaNumber(4);
//  auto volumeFour=ROhist-> GetVolume(4)->GetName();
//
//  auto stageFive=ROhist->GetReplicaNumber(5);
//  auto volumeFive=ROhist-> GetVolume(5)->GetName();

  //G4cout<<stageNull<<volumeNull<<"\t"<<stageOne<<volumeOne<<"\t"<<stageTwo<<volumeTwo<<"\t"<<stageThree<<volumeThree<<"\t"<</*stageFour<<stageFive<<*/G4endl;

 }


  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[layerNumber];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber; 
    G4Exception("B4cCalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }         

  // Get hit for total accounting
  auto hitTotal 
    = (*fHitsCollection)[fHitsCollection->entries()-1];
  
  // Add values
  hit->Add(edep, stepLength);
  hitTotal->Add(edep, stepLength); 
      
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
