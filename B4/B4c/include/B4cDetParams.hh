#ifndef B4cDetParams
#define B4cDetParams

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"



//Specifyes the geometrical information

class DetParams {

public:

static DetParams& Instance(){

								static DetParams Det;
								return Det;

}



~DetParams();

void InitDet();

void SetfNofLayers(G4double nla);
G4double GetfNofLayers();

void SettileLenX(G4double tx);
G4double GettileLenX();

void SettileLenY(G4double ty);
G4double GettileLenY();

void SetabsoThickness(G4double abs);
G4double GetabsoThickness();

void SetgapThickness(G4double gap);
G4double GetgapThickness();

void SetcalorSizeXY(G4double cs);
G4double GetcalorSizeXY();

void SetWorldMult(G4double wm);

G4double GetnofTilesX();
G4double GetnofTilesY();
G4double GetlayerThickness();
G4double GetcalorThickness();
G4double GettilesPerLayer();
G4double GetWorldSizeXY();
G4double GetWorldSizeZ();
//G4double GetGunPos();

private:

G4double fNofLayers;

G4double tileLenX;
G4double tileLenY;
G4double nofTilesX;
G4double nofTilesY;
G4double tilesPerLayer;


G4double absoThickness;
G4double gapThickness;
G4double layerThickness;
G4double calorSizeXY;
G4double calorThickness;

G4double WorldSizeXY;
G4double WorldSizeZ;
G4double WorldMult;

G4double GunPos;



DetParams(){
}
DetParams(const DetParams& other){
}

};


inline DetParams& GetInst(){
								return DetParams::Instance();
}

#endif
