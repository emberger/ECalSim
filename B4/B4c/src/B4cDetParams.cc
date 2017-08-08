#include "B4cDetParams.hh"



DetParams::~DetParams(){
}

void DetParams::InitDet(){

								nofTilesX=calorSizeXY/tileLenX;
								nofTilesY=calorSizeXY/tileLenY;

								layerThickness=absoThickness+gapThickness;

								calorThickness=fNofLayers*layerThickness;

								tilesPerLayer=nofTilesX*nofTilesY;

								WorldSizeXY=2.0 * calorSizeXY;
								WorldSizeZ=WorldMult * calorThickness;
}

void DetParams::SetfNofLayers(G4double nla){
								fNofLayers=nla;
}

G4double DetParams::GetfNofLayers(){
								return fNofLayers;
}

void DetParams::SettileLenX(G4double tx){
								tileLenX=tx;
}

G4double DetParams::GettileLenX(){
								return tileLenX;
}

void DetParams::SettileLenY(G4double ty){
								tileLenY=ty;
}

G4double DetParams::GettileLenY(){
								return tileLenY;
}

void DetParams::SetabsoThickness(G4double abso){
								absoThickness=abso;
}

G4double DetParams::GetabsoThickness(){
								return absoThickness;
}

void DetParams::SetgapThickness(G4double gap){
								gapThickness=gap;
}

G4double DetParams::GetgapThickness(){
								return gapThickness;
}

void DetParams::SetcalorSizeXY(G4double cs){
								calorSizeXY=cs;
}

G4double DetParams::GetcalorSizeXY(){
								return calorSizeXY;
}

void DetParams::SetWorldMult(G4double wm){
								WorldMult=wm;
}



G4double DetParams::GetnofTilesX(){
								return nofTilesX;
}
G4double DetParams::GetnofTilesY(){
								return nofTilesY;
}
G4double DetParams::GetlayerThickness(){
								return layerThickness;
}
G4double DetParams::GetcalorThickness(){
								return calorThickness;
}
G4double DetParams::GettilesPerLayer(){
								return tilesPerLayer;
}
G4double DetParams::GetWorldSizeXY(){
								return WorldSizeXY;
}
G4double DetParams::GetWorldSizeZ(){
								return WorldSizeZ;
}
