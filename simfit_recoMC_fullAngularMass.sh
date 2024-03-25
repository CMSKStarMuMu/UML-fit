#!/bin/bash

par=1

## full stat
multi=0
nsam=0
#### data like one job
# multi=0
# nsam=${1}

plot=1
save=2

XGBv=8

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit4d ]; then mkdir logs_simFit4d; fi
if [ ! -d simFitResults4d ]; then mkdir simFitResults4d; fi
if [ ! -d plotSimFit4d_d ]; then mkdir plotSimFit4d_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullAngularMass; then

    while read -a line; do
	bin=${line[0]}

	for year in {2016..2018}; do

	    ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} ${XGBv} 0 ${plot} ${save} ${year} \
		&>logs_simFit4d/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_${XGBv}_${year}.out &

	done

	./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} ${XGBv} 0 ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit4d/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_${XGBv}_2016_2017_2018.out &

    done < ../confSF/KDE_SF.list

fi
