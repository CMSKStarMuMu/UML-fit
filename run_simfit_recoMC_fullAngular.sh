#!/bin/bash

par=1

multi=0
nsam=${1}

plot=0
save=1

ibin=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-custMinos
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export EFFDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/files

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
if [ "${par}" == 0 ]; then
    parstr="ev"
else
    parstr="od"
fi
if [ ! -r $EFFDIR/KDEeff_b${bin}_${parstr}_2016.root ]; then
    echo $EFFDIR/KDEeff_b${bin}_${parstr}_2016.root not found
    exit 1
fi
if [ ! -r $EFFDIR/KDEeff_b${bin}_${parstr}_2017.root ]; then
    echo $EFFDIR/KDEeff_b${bin}_${parstr}_2017.root not found
    exit 1
fi
if [ ! -r $EFFDIR/KDEeff_b${bin}_${parstr}_2018.root ]; then
    echo $EFFDIR/KDEeff_b${bin}_${parstr}_2018.root not found
    exit 1
fi
if [ ! -r $HOME/simfit_recoMC_fullAngular ]; then
    echo $HOME/simfit_recoMC_fullAngular not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root .
cp $EFFDIR/KDEeff_b${bin}_${parstr}_2016.root .
cp $EFFDIR/KDEeff_b${bin}_${parstr}_2017.root .
cp $EFFDIR/KDEeff_b${bin}_${parstr}_2018.root .
cp $HOME/simfit_recoMC_fullAngular .

mkdir simFitResults
mkdir plotSimFit_d

echo ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018
./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} 1 ${plot} ${save} 2016 2017 2018

if [ ! -d $HOME/simFitResults ]; then
    mkdir $HOME/simFitResults
fi
if [ ! -d $HOME/plotSimFit_d ]; then
    mkdir $HOME/plotSimFit_d
fi
cp plotSimFit_d/* $HOME/plotSimFit_d/
cp simFitResults/* $HOME/simFitResults/
# for file in simFitResults/* ; do cp $file $HOME/${file//.root/_${multi}s${nsam}.root}; done

rm -rf plotSimFit_d
rm -rf simFitResults

rm simfit_recoMC_fullAngular
rm recoMCDataset_b*
rm KDEeff_b*
rm simFitResult*
