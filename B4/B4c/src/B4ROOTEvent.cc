// Implementation file for Pi0Event and components
// create dictionary with rootcint RootClasses_dict.cpp -c include/B4ROOTEvent.hh include/LinkDef.h
//

#include "B4ROOTEvent.hh"
//#include "Stiostream.h"

ClassImp(B4ROOTEvent)
ClassImp(B4ROOTHit)


TClonesArray* B4ROOTEvent::aHits = 0;

// TPi0Event class

B4ROOTEvent::B4ROOTEvent() {
	if (!aHits) aHits = new TClonesArray("B4ROOTHit", 1000);
	m_Hits = aHits;
	m_EventNo = 0;
	m_GapEnergy = 0;
	m_NHits = 0;

}

B4ROOTEvent::~B4ROOTEvent() {
	Clear();
}

void B4ROOTEvent::Clear(const Option_t*) {
	std::cout << "Clearing TPi0Event..." << std::endl;
	m_EventNo = 0;
	m_GapEnergy = 0;
	m_NHits=0;
	m_Hits->Clear();
	std::cout << "Clear done!" << std::endl;
}
void B4ROOTEvent::Reset(Option_t * /*option*/)
{
	// Static function to reset all static objects for this event
		delete aHits; aHits = 0;
}




B4ROOTHit* B4ROOTEvent::AddHit(B4ROOTHit& cand) {
	//cout << "Adding Pi0 candidate after " << m_NPi0Candidates << endl;
	TClonesArray &aCand = *m_Hits;
	B4ROOTHit *piCand = new(aCand[m_NHits++]) B4ROOTHit(cand);
	//cout << "Added Candidate " << m_NPi0Candidates << endl;
	return piCand;
}

B4ROOTHit::B4ROOTHit() : TObject() {
	m_X = 0;
	m_Y = 0;
	m_Z = 0;
	m_EnergyDeposit = 0;
}

B4ROOTHit::B4ROOTHit(const B4ROOTHit& orig) : TObject(orig) {
	m_X = orig.m_X;
	m_Y = orig.m_Y;
	m_Z = orig.m_Z;
	m_EnergyDeposit = orig.m_EnergyDeposit;
}

B4ROOTHit::~B4ROOTHit(){}
