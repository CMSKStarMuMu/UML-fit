#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
Executable  = run_simfit_data_fullAngularMass_Swave_toyeff.sh
toy	    = \$(ProcId) % 20
bin         = ( \$(ProcId) / 20 )
Arguments   = \$INT(toy) \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_data_fullAngularMass_Swave_toyeff_\$INT(toy)_\$INT(bin).out
Error       = logs_parSub/simfit_data_fullAngularMass_Swave_toyeff_\$INT(toy)_\$INT(bin).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
echo 'Queue 120'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# Compilation, submission and file removal
if make simfit_data_fullAngularMass_Swave_toyeff
then condor_submit temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
rm temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# yearConf    = \$(ProcId) % 4
