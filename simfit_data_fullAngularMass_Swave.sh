#!/bin/bash

par=1

multi=0
nsam=0
q2stat=0

localFile=0

plot=1
save=2

bin=${1}
yearConf=${2}
XGBv=${3}
[ -z "${XGBv}" ] && XGBv=0	# set default
fitopt=${4}
[ -z "${fitopt}" ] && fitopt=0	# set default

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit4d ]; then mkdir logs_simFit4d; fi
if [ ! -d simFitResults4d ]; then mkdir simFitResults4d; fi
if [ ! -d plotSimFit4d_d ]; then mkdir plotSimFit4d_d; fi


# Compile dictionary and macro
# make AngDict
if make simfit_data_fullAngularMass_Swave; then

    case "$yearConf" in

	0)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2016 2017 2018
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2016 2017 2018 \
		&> logs_simFit4d/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_${XGBv}_${fitopt}_${localFile}_${plot}_${save}_2016_2017_2018.log &
	    ;;

	1)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2016
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2016 \
		&> logs_simFit4d/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_${XGBv}_${fitopt}_${localFile}_${plot}_${save}_2016.log &
	    ;;

	2)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2017
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2017 \
		&> logs_simFit4d/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_${XGBv}_${fitopt}_${localFile}_${plot}_${save}_2017.log &
	    ;;

	3)
	    echo ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2018
	    ./simfit_data_fullAngularMass_Swave ${bin} ${par} ${multi} ${nsam} ${q2stat} ${XGBv} ${fitopt} ${localFile} ${plot} ${save} 2018 \
		&> logs_simFit4d/simfit_data_fullAngularMass_Swave_${bin}_${par}_${multi}_${nsam}_${q2stat}_${XGBv}_${fitopt}_${localFile}_${plot}_${save}_2018.log &
	    ;;

    esac

fi
