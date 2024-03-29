#!/bin/bash

seed=${1}
nsam=10

save=1

bin=${2}
gen1=${3}
gen2=${4}
gen3=${5}
gen4=${6}
gen5=${7}
gen6=${8}
gen7=${9}
gen8=${10}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-custMinos
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

if [ ! -d $HOME/toyFitResults_b${bin} ]; then
    mkdir $HOME/toyFitResults_b${bin}
fi
if [ ! -d $HOME/plotSimFit_b${bin} ]; then
    mkdir $HOME/plotSimFit_b${bin}
fi

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
if [ ! -r $HOME/simfit_toy_fullAngular ]; then
    echo $HOME/simfit_toy_fullAngular not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root .
cp $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root .
cp $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root .
cp $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root .
cp $HOME/simfit_toy_fullAngular .

mkdir toyFitResults
mkdir plotSimFit_d

echo ./simfit_toy_fullAngular ${bin} ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8} ${seed} ${nsam} 1 ${save} 2016 2017 2018
./simfit_toy_fullAngular ${bin} ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8} ${seed} ${nsam} 1 ${save} 2016 2017 2018

cp toyFitResults/* $HOME/toyFitResults_b${bin}/
cp plotSimFit_d/* $HOME/plotSimFit_b${bin}/

rm -rf toyFitResults
rm -rf plotSimFit_d

rm simfit_toy_fullAngular
rm recoMCDataset_b*
rm KDEeff_b*
