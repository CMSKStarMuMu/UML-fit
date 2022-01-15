#!/bin/bash

par=0
sample=${1}

pow=${2}
co1=${3}
co4=${4}
co5=${5}

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
make simfit_genMC

while read -a line; do
    bin=${line[0]}

    nohup ./simfit_genMC ${bin} ${par} ${sample} ${pow} ${co1} ${co4} ${co5} 1 1 \
	&> logs_simFit/simfit_genMC_${bin}_${par}_${sample}_${pow}_${co1}_${co4}_${co5}.out &

done < ../confSF/KDE_SF.list
