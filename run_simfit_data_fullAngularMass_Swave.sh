#!/bin/bash

par=1

multi=0
nsam=${1}
q2stat=${4}

XGBv=8
localFile=0
fitopt=0
unbl=4

plot=1
save=${5}

bin=${2}
# ibin=${2}

yearConf=${3}

export SAMPLEDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/data-datasets/
export SBDIR=/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/sidebands/

if [ "${XGBv}" == 0 ]; then
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
        dataname="${SAMPLEDIR}/recoDATADataset_b${bin}_${iy}.root"
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

if [ ! -r $HOME/simfit_data_fullAngularMass_Swave ]; then
    echo $HOME/simfit_data_fullAngularMass_Swave not found
    exit 1
fi
cp $HOME/simfit_data_fullAngularMass_Swave .
cp $HOME/*.pcm .

mkdir simFitResults4d
mkdir plotSimFit4d_d

case "$yearConf" in
# void simfit_data_fullAngularMass_SwaveBin(int q2Bin, int parity, bool multiSample, uint nSample, uint q2stat, int fitOption, int XGBv, int unblind, bool localFiles, bool plot, int save, std::vector<int> years)

    0)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016 2017 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016 2017 2018
	;;

    1)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2016
	;;

    2)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2017
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2017
	;;

    3)
	echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2018
	./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${fitopt} ${XGBv} ${unbl} ${localFile} ${plot} ${save} 2018
  ;;

esac

if [ ! -d $HOME/simFitResults4d ]; then
    mkdir $HOME/simFitResults4d
fi
if [ ! -d $HOME/plotSimFit4d_d ]; then
    mkdir $HOME/plotSimFit4d_d
fi
cp plotSimFit4d_d/* $HOME/plotSimFit4d_d/
cp simFitResults4d/* $HOME/simFitResults4d/
