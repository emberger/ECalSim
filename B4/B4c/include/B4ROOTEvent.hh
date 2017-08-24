#ifndef B4ROOTEvent_hh
#define B4ROOTEvent_hh

// Header file for Pi0 Event class and sub-classes with specific information
// create dictionary with rootcint RootClasses_dict.cpp -c include/B4ROOTEvent.hh include/LinkDef.h


#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TMath.h"
#include "TVector3.h"
#include <iostream>

//#include "globals.hh"

class B4ROOTHit : public TObject {

private:
Double_t m_X;
Double_t m_Y;
Double_t m_Z;
Double_t m_EnergyDeposit;

Int_t m_PhotonNr;

public:
B4ROOTHit();
B4ROOTHit(const B4ROOTHit& orig);
virtual ~B4ROOTHit();
void Clear(const Option_t* option ="") {
};
// Getters:
Double_t X() const {
								return m_X;
}
Double_t Y() const {
								return m_Y;
}
Double_t Z() const {
								return m_Z;
}
Double_t EnergyDeposit() const {
								return m_EnergyDeposit;
}

Int_t PhotNr(){
								return m_PhotonNr;
}

// Setters:
void SetPhotNr(Int_t nr){
								m_PhotonNr=nr;
}

void SetCoordinates(Double_t x, Double_t y, Double_t z) {
								m_X = x; m_Y = y; m_Z = z;
}
void SetEnergyDeposit(Double_t dep) {
								m_EnergyDeposit = dep;
}

ClassDef(B4ROOTHit,1)
};




class B4ROOTEvent : public TObject {

private:
Int_t m_EventNo;

Double_t m_GapEnergy;

Int_t m_NHits;

Double_t m_AbsoThickness;
Double_t m_GapThickness;
Double_t m_Layerno;
Double_t m_tilesizeX;
Double_t m_tilesizeY;
Double_t m_calsizeXY;

Double_t m_EnergyPhoton1;
Double_t m_EnergyPhoton2;

TVector3 m_MomentumPh1;
TVector3 m_MomentumPh2;


TClonesArray *m_Hits;  //->
static TClonesArray *aHits;

public:
B4ROOTEvent();
virtual ~B4ROOTEvent();
void Clear(const Option_t* option ="");
static void   Reset(Option_t *option ="");



//Getters
Int_t EventNo() const {
								return m_EventNo;
}
Double_t GapEnergy() const {
								return m_GapEnergy;
}

Double_t AbsoThickness(){
								return m_AbsoThickness;
}
Double_t GapThickness(){
								return m_GapThickness;
}
Double_t NumberOfLayers(){
								return m_Layerno;
}
Double_t TilesizeX(){
								return m_tilesizeX;
}
Double_t TilesizeY(){
								return m_tilesizeY;
}
Double_t calsizeXY(){
								return m_calsizeXY;
}

Double_t EnergyPhoton1(){
								return m_EnergyPhoton1;
}
Double_t EnergyPhoton2(){
								return m_EnergyPhoton2;
}

TVector3 MomentumPh1(){
								return m_MomentumPh1;
}

TVector3 MomentumPh2(){
								return m_MomentumPh2;
}


Int_t NHits() const {
								return m_NHits;
}

TClonesArray* Hits() const {
								return m_Hits;
}

B4ROOTHit* Hit(Int_t i) {
								return (B4ROOTHit*)m_Hits->At(i);
}

//Setters
void SetEventNo(Int_t evN) {
								m_EventNo = evN;
}
void SetGapEnergy(Double_t en) {
								m_GapEnergy = en;
}

void SetAbsoThickness(Double_t at){
								m_AbsoThickness=at;
}
void SetGapThickness(Double_t gt){
								m_GapThickness=gt;
}
void SetNumberOfLayers(Double_t nol){
								m_Layerno=nol;
}
void SetTilesizeX(Double_t tsX){
								m_tilesizeX=tsX;
}
void SetTilesizeY(Double_t tsY){
								m_tilesizeY=tsY;
}
void SetcalsizeXY(Double_t cxy){
								m_calsizeXY=cxy;
}

void SetEnergyPhoton1(Double_t eph1){
								m_EnergyPhoton1=eph1;
}

void SetEnergyPhoton2(Double_t eph2){
								m_EnergyPhoton2=eph2;
}

void SetMomentumPh1(Double_t mx1, Double_t my1, Double_t mz1){
								m_MomentumPh1.SetXYZ(mx1, my1, mz1);

}

void SetMomentumPh2(Double_t mx2, Double_t my2, Double_t mz2){
								m_MomentumPh2.SetXYZ(mx2, my2, mz2);
}

B4ROOTHit* AddHit(B4ROOTHit& cand);

ClassDef(B4ROOTEvent,1)
};
#endif
