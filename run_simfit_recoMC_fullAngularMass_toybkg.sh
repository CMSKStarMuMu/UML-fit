#!/bin/bash
multi=0
xgb=8

plot=1
save=1

nsam=${1}
ibin=${2}
par=${3}

## if localFile == 1, please check the name of the MC dataset being copied
localFile=0

export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-datasets/
export SBDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/

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
        echo 'will copy data samples from  from ' ${SAMPLEDIR}
        echo 'will copy efficiencies from ' ${EFFDIR}
        echo 'will copy sidebands from ' ${SBDIR}
        [ "$yearConf" -gt 0 ] && [ "$((${yearConf}+2015))" != "$iy" ] && continue
        if [ "${xgb}" == 0 ]; then
            ## this does not exists at the moment
            dataname="${SAMPLEDIR}/recoMCDataset_b${bin}_${iy}.root"  
        else    
            dataname="${SAMPLEDIR}/recoMCDataset_b${bin}_${iy}_XGBv8.root"
        fi    
        effname="${EFFDIR}/KDEeff_b${bin}_${parstr}_${iy}.root"
        sbname="${SBDIR}/savesb_${iy}_b${bin}.root"
        if [ ! -r "${dataname}" ]; then
            echo "${dataname}" not found
            exit 1
        fi
        if [ ! -r "${effname}" ]; then
            echo "${effname}" not found
            exit 1
        fi
        if [ ! -r "${sbname}" ]; then
            echo "${sbname}" not found
            exit 1
        fi
        cp "${dataname}" .
        cp "${effname}" .
        cp "${sbname}" .
    done
fi

if [ ! -r $HOME/simfit_recoMC_fullAngularMass_toybkg ]; then
    echo $HOME/simfit_recoMC_fullAngularMass_toybkg not found
    exit 1
fi
cp $HOME/simfit_recoMC_fullAngularMass_toybkg .
cp $HOME/*.pcm .

mkdir -p simFitResults4d/xgbv8
mkdir -p plotSimFit4d_d/xgbv8

if [ "$bin" -lt 8 ]; then
    echo ./simfit_recoMC_fullAngularMass_toybkg ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 2017 2018
    ./simfit_recoMC_fullAngularMass_toybkg ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 2017 2018
else
    echo ./simfit_recoMC_fullAngularMass_toybkg ${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016
    ./simfit_recoMC_fullAngularMass_toybkg${bin} ${par} ${multi} ${nsam} ${xgb} ${localFile} ${plot} ${save} 2016 
fi


if [ ! -d $HOME/simFitResults4d/xgbv8/ ]; then
    mkdir -p $HOME/simFitResults4d/xgbv8/
fi
if [ ! -d $HOME/plotSimFit4d_d/xgbv8/ ]; then
    mkdir -p $HOME/plotSimFit4d_d/xgbv8/
fi
cp plotSimFit4d_d/xgbv8/* $HOME/plotSimFit4d_d/xgbv8/
cp simFitResults4d/xgbv8/* $HOME/simFitResults4d/xgbv8/
