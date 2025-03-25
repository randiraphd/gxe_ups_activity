
#Shell script to prep sequencing files on the cluster

## -----
## dir w/ raw reads
raw=/yoursequencingdirectory/
ls -alh ${raw} #shows files and details in the folder named 'raw'


## -----
## make project sub-directories
project=/myproteasome/data/illumina/2023.04_GxE/
mkdir -pv ~/${project}{scripts,alignments,logs,error_logs,std_err_log}

#now add map_filter_count script to the scripts folder

ls -alh ~/${project}
