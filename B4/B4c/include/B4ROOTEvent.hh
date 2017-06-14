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
#include <iostream>

//#include "globals.hh"

class B4ROOTHit : public TObject {

private:
	Double_t m_X;
	Double_t m_Y;
	Double_t m_Z;
	Double_t m_EnergyDeposit;


public:
	B4ROOTHit();
	B4ROOTHit(const B4ROOTHit& orig);
	virtual ~B4ROOTHit();
	void Clear(const Option_t* option ="") {};
	// Getters:
	Double_t X() const {return m_X;}
	Double_t Y() const {return m_Y;}
	Double_t Z() const {return m_Z;}
	Double_t EnergyDeposit() const {return m_EnergyDeposit;}

	// Setters:
	void SetCoordinates(Double_t x, Double_t y, Double_t z) {m_X = x; m_Y = y; m_Z = z;}
	void SetEnergyDeposit(Double_t dep) {m_EnergyDeposit = dep;}

	ClassDef(B4ROOTHit,3)
};




class B4ROOTEvent : public TObject {

private:
	Int_t m_EventNo;
	Double_t m_GapEnergy;
	Int_t m_NHits;
	Int_t m_TilesX;
	Int_t m_TilesY;
	Int_t m_Layers;
	TClonesArray *m_Hits; //->
	static TClonesArray *aHits;

public:
	B4ROOTEvent();
	virtual ~B4ROOTEvent();
	void Clear(const Option_t* option ="");
	static void   Reset(Option_t *option ="");



	//Getters
	Int_t EventNo() const {return m_EventNo;}
	Double_t GapEnergy() const {return m_GapEnergy;}
	Int_t TilesX(){return m_TilesX;}
	Int_t TilesY(){return m_TilesY;}
	Int_t Layers(){return m_Layers;}

	Int_t NHits() const {return m_NHits;}
	TClonesArray* Hits() const {return m_Hits;}
	B4ROOTHit* Hit(Int_t i) {return (B4ROOTHit*)m_Hits->At(i);}

	//Setters
	void SetEventNo(Int_t evN) {m_EventNo = evN;}
	void SetGapEnergy(Double_t en) {m_GapEnergy = en;}
	void SetTilesX(Int_t nx){m_TilesX=nx;}
	void SetTilesY(Int_t ny){m_TilesY=ny;}
	void SetLayers(Int_t nl){m_Layers=nl;}


	B4ROOTHit* AddHit(B4ROOTHit& cand);

	ClassDef(B4ROOTEvent,1)
};
#endif
