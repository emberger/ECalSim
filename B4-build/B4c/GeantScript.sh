#!/bin/bash

echo working!


AngCOUNTER=0


Sen=100          #starting energy
EnIterations=5
increment=100
eventsPerEnergy=1000

Xangle=-1.0      #starting angles
Yangle=-1.0
Zangle=1.0
XAngIncrement=0.2
YAngIncrement=0.2
AngIiterations=10





foldername=MuonEnergyandSlopeScan_5x5mmTile_1mmLead   #folder containing .root files






mkdir $foldername

sed -i "19 s%^/run/beamOn.*%/run/beamOn $eventsPerEnergy%" GeantSim.mac

    while [ $AngCOUNTER -lt $AngIiterations ]; do

        sed -i "16 s%^/gun/direction.*%/gun/direction $Xangle $Yangle $Zangle%" GeantSim.mac

        EnCOUNTER=0
        en=$Sen
        buf=start
        sed -i "18 s%^/gun/energy.*%/gun/energy $buf MeV%" GeantSim.mac

        while [  $EnCOUNTER -lt $EnIterations ]; do

            s=$foldername

            sed -i "18 s%$buf%$en%" GeantSim.mac
            buf=$en

            ./exampleB4c -m GeantSim.mac

            s+=/Muon
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
