#!/bin/bash
echo $SGE_TASK_ID
mkdir ../$SGE_TASK_ID
cd ../$SGE_TASK_ID
cp ../comm/exe ./
cp ../comm/insulator_BN_test.py ./
cp ../data/ELoss_hBN.txt ./
cp ../data/DESCS_BN.txt ./
cp ../data/ThetaEl_BN.txt ./
./exe
cd ../comm
