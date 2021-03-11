#!/bin/bash

par=1

multi=0
nsam=${1}

plot=1
save=1

ibin=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-integration
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export SAMPLEDIR2=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

nbin=0
while read -a line; do
    abin[$nbin]=${line[0]}
    nbin=$((nbin+1))
done < $HOME/../confSF/KDE_SF.list
bin=${abin[$ibin]}

echo 'now submitting for bin ' ${bin}

if [ ! -r $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR2/2016/lmnr/integrSwave_b${bin}_od_2016.root ]; then
    echo $SAMPLEDIR2/2016/lmnr/integrSwave_b${bin}_od_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR2/2017/lmnr/integrSwave_b${bin}_od_2017.root ]; then
    echo $SAMPLEDIR2/2017/lmnr/integrSwave_b${bin}_od_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR2/2018/lmnr/integrSwave_b${bin}_od_2018.root ]; then
    echo $SAMPLEDIR2/2018/lmnr/integrSwave_b${bin}_od_2018.root not found
    exit 1
fi
if [ ! -r $HOME/simfit_recoMC_fullAngular_Swave ]; then
    echo $HOME/simfit_recoMC_fullAngular_Swave not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root .
cp $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root .
cp $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root .
cp $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root .
cp $SAMPLEDIR2/2016/lmnr/integrSwave_b${bin}_od_2016.root .
cp $SAMPLEDIR2/2017/lmnr/integrSwave_b${bin}_od_2017.root .
cp $SAMPLEDIR2/2018/lmnr/integrSwave_b${bin}_od_2018.root .
cp $HOME/simfit_recoMC_fullAngular_Swave .

mkdir simFitResults
mkdir plotSimFit_d

echo ./simfit_recoMC_fullAngular_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
./simfit_recoMC_fullAngular_Swave ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018

if [ ! -d $HOME/simFitResults ]; then
    mkdir $HOME/simFitResults
fi
if [ ! -d $HOME/plotSimFit_d ]; then
    mkdir $HOME/plotSimFit_d
fi
cp plotSimFit_d/* $HOME/plotSimFit_d/
cp simFitResults/* $HOME/simFitResults/

rm -rf plotSimFit_d
rm -rf simFitResults

rm simfit_recoMC_fullAngular_Swave
rm recoMCDataset_b*
rm KDEeff_b*
rm integrSwave_b*
rm simFitResult*
