#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir -p logs_parSub; fi

nbins=$(wc -l ../confSF/KDE_SF_all.list | cut -d " " -f1)

while read -a line; do
    bin=${line[0]}
    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_simfit_recoMC_fullAngular_oneBin${bin}.sub
Executable  = run_simfit_recoMC_fullAngular.sh
nsamp       = \$(ProcId)
parity      = 1
yearConf    = 0
Arguments   = \$INT(nsamp) ${bin} \$INT(parity) \$INT(yearConf)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_recoMC_fullAngular_\$INT(nsamp)_p\$INT(parity)_b${bin}.out
Error       = logs_parSub/simfit_recoMC_fullAngular_\$INT(nsamp)_p\$INT(parity)_b${bin}.err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

    if [ "${USER}" == "fiorendi" ]; then
        echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngular_oneBin{bin}.sub
    fi
    echo "Queue 100">>temp_sub_simfit_recoMC_fullAngular_oneBin${bin}.sub

    # Compilation, submission and file removal
    if make simfit_recoMC_fullAngular
    then condor_submit temp_sub_simfit_recoMC_fullAngular_oneBin${bin}.sub
    fi
    rm temp_sub_simfit_recoMC_fullAngular_oneBin{bin}.sub
done < ../confSF/KDE_SF_all.list
