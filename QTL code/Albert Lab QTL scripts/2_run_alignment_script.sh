## -----
## log in to cluster


## -----
## setup
dir=/home/albertf/randia/myproteasome/data/illumina/
proj=${dir}2023.04_GxE/
script=${proj}/scripts/map_filter_count_4cluster.sh
tail ${script}
output_log=${proj}logs/std_out_log
err_log=${proj}logs/std_err_log


## -----
## make script executable:
chmod 777 ${script}


## -----
## run the script via 'sbatch' and direct output to log:
## '&' runs the script in the background...
sbatch ${script} & >> ${output_log}


## -----
## read logs and check output:
cat ${output_log} #bash shell utilities
cat ${err_log}
ls -alh ${proj}/alignments

## make sure your script is still running:
#squeue -u x500 tells you if you have a job running
