#!/bin/bash
multi=0
xgb=8

plot=0
save=1

nsam=${1}
ibin=${2}
par=${3}
yearConf=${4}

## if localFile == 1, please check the name of the MC dataset being copied
localFile=0

export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-datasets/

if [ "${xgb}" == 0 ]; then
    export EFFDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-theta-v7/files
else
    export EFFDIR=/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-theta-v7-XGBv8/files
fi


if [ "${USER}" == "aboletti" ]; then
    export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-JpsiFit
    export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
elif [ "${USER}" == "fiorendi" ]; then
    export HOME=/afs/cern.ch/work/f/fiorendi/private/effKDE/UML-fit
    export CMSSWDIR=/afs/cern.ch/work/f/fiorendi/private/effKDE/CMSSW_10_4_0/src
else
    echo no user found
    exit 1
fi

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

# nbin=0
# while read -a line; do
#     abin[$nbin]=${line[0]}
#     nbin=$((nbin+1))
# done < $HOME/../confSF/KDE_SF.list
# bin=${abin[$ibin]}
bin=${ibin}

echo 'now submitting for bin ' ${bin}

if [ "${par}" == 0 ]; then
    parstr="ev"
else
    parstr="od"
fi

if [ "${localFile}" -gt 0 ]; then 
    for iy in {2016..2018}
    do
        echo 'will copy MC samples from  from ' ${SAMPLEDIR}
        echo 'will copy efficiencies from ' ${EFFDIR}
        if [ "${xgb}" == 0 ]; then
            ## this does not exists at the moment
            dataname="${SAMPLEDIR}/recoMCDataset_b${bin}_${iy}.root"  
        else    
            dataname="${SAMPLEDIR}/recoMCDataset_b${bin}_${iy}_XGBv8.root"
        fi    
        effname="${EFFDIR}/KDEeff_b${bin}_${parstr}_${iy}.root"
        if [ ! -r "${dataname}" ]; then
            echo "${dataname}" not found
            exit 1
        fi
        if [ ! -r "${effname}" ]; then
            echo "${effname}" not found
            exit 1
        fi
        cp "${dataname}" .
        cp "${effname}" .
    done
fi


if [ ! -r $HOME/simfit_recoMC_fullAngular ]; then
    echo $HOME/simfit_recoMC_fullAngular not found
    exit 1
fi
cp $HOME/simfit_recoMC_fullAngular .

mkdir -p simFitResults/xgbv8
mkdir -p plotSimFit_d/xgbv8

case "$yearConf" in
    0)
        echo ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 2017 2018
        ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 2017 2018
	;;

    1)
        echo ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016
        ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 
	;;

    2)
        echo ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2017
        ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2017 
	;;

    3)
        echo ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2018
        ./simfit_recoMC_fullAngular ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2018 
  ;;

esac

if [ ! -d $HOME/simFitResults/xgbv8 ]; then
    mkdir -p $HOME/simFitResults/xgbv8
fi
if [ ! -d $HOME/plotSimFit_d/xgbv8 ]; then
    mkdir -p $HOME/plotSimFit_d/xgbv8
fi
cp plotSimFit_d/xgbv8/* $HOME/plotSimFit_d/xgbv8/
cp simFitResults/xgbv8/* $HOME/simFitResults/xgbv8/

# for file in simFitResults/* ; do cp $file $HOME/${file//.root/_${multi}s${nsam}.root}; done

# rm -rf plotSimFit_d
# rm -rf simFitResults
# 
# rm simfit_recoMC_fullAngular
# rm recoMCDataset_b*
# rm KDEeff_b*
# rm simFitResult*
