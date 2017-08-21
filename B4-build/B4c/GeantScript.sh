#!/bin/bash

echo working!


AngCOUNTER=0


Sen=100          #starting energy
EnIterations=10
increment=100
eventsPerEnergy=2000

Xangle=-0.2      #starting angles
Yangle=-0.2
Zangle=1.0
XAngIncrement=-0.2
YAngIncrement=-0.2
AngIiterations=5





foldername=GammaEnergyandSlopeScan   #folder containing .root files






mkdir $foldername

sed -i "19 s%^/run/beamOn.*%/run/beamOn $eventsPerEnergy%" run2.mac

    while [ $AngCOUNTER -lt $AngIiterations ]; do

        sed -i "16 s%^/gun/direction.*%/gun/direction $Xangle $Yangle $Zangle%" run2.mac

        EnCOUNTER=0
        en=$Sen
        buf=start
        sed -i "18 s%^/gun/energy.*%/gun/energy $buf MeV%" run2.mac

        while [  $EnCOUNTER -lt $EnIterations ]; do

            s=$foldername

            sed -i "18 s%$buf%$en%" run2.mac
            buf=$en

            ./exampleB4c -m GeantSim.mac

            s+=/Gamma
            s+=$en
            s+=MeV
            s+=_$Xangle
            s+=:$Yangle
            s+=_.root

            #echo $s

            mv Tree.root ./$s

            let en=en+$increment
            let EnCOUNTER=EnCOUNTER+1
        done


        Xangle=$(echo "($Xangle+$XAngIncrement)" | bc)
        Yangle=$(echo "($Yangle+$YAngIncrement)" | bc)

        let AngCOUNTER=AngCOUNTER+1
done
