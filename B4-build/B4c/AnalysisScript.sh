
#!/bin/bash

Path=/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/

FILES=/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/GammaEnergyandSlopeScan_5x5mmTile_1mmLead/*.root

mkdir $Path/GammaEnergyandSlopeScan_5x5mmTile_1mmLead_analysis

AnaPath=$Path/GammaEnergyandSlopeScan_5x5mmTile_1mmLead_analysis


counter=0
for f in $FILES
do
    #echo "$f"
    FilePath=$f
    filename="${FilePath##*/}"

    foldername=${filename%_*}

    mkdir $AnaPath/$foldername

    ./Analysis 1 50 0 1000 0 $FilePath $AnaPath/$foldername


done
