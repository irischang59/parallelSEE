#!/bin/bash
for i in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99 
do
mkdir $i
cp -r comm ./$i
cd $i/comm
sed -i '3s/.*/python insulator_BN.py '$i' 1.75e-06 100 >> hBNSEY.txt/' exe
qsub -cwd -V -N hBN_0_$i -l h_data=2G,h_rt=8:00:00 -t 1-10:1 /u/home/i/irischan/hBN_0/$i/comm/run.sh
cd ../../
done
