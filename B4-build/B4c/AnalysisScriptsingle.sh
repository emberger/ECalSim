
#!/bin/bash

Path=/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c

FILES=/home/iwsatlas1/emberger/Geant4/Current/SensitiveDetector/B4-build/B4c/rmsx_over_E/10cm/02/gamma300MeV_0.2:0.2_.root

#mkdir $Path/rmsx_over_E_analysis/10cm/02

#AnaPath=$Path/rmsx_over_E_analysis/10cm/02


counter=0
for f in $FILES
do
    #echo "$f"
    FilePath=$f
    filename="${FilePath##*/}"

    foldername=${filename%_*}

    mkdir $AnaPath/$foldername

    ./Analysis 1 50 0 1000 0 $FilePath $AnaPath/$foldername 102


done
