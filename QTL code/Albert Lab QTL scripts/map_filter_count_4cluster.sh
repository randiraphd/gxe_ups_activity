#!/bin/bash -l
#SBATCH --time=35:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p ram256g
#SBATCH -n 2
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=randia@umn.edu
#SBATCH -o /home/albertf/randia/myproteasome/data/illumina/2023.04_GxE/logs/std_out_log
#SBATCH -e /home/albertf/randia/myproteasome/data/illumina/2023.04_GxE/logs/std_err_log


## -----
## load software; the '/xx' denotes a version
bwaVersion='bwa/0.7.15'
module load ${bwaVersion}
samtoolsVersion='samtools/1.5'
module load ${samtoolsVersion}
module load bamtools


## -----

## specify raw data and output location:
raw=/yoursequencingdirectory/

ls -alh ${raw}

## specify output location (logs, alignments, etc....):
out=/home/albertf/randia/myproteasome/data/illumina/2023.04_GxE/
ls -alh ${out}

BWA_data=${out}alignments
ls -alh ${BWA_data}

output_log=/home/albertf/randia/myproteasome/data/illumina/2023.04_GxE/logs/std_out_log
touch ${output_log}
echo "" > ${output_log}
echo -e "On $( date ) started alignments:\n" >> ${output_log}
cat ${output_log}


## make sure we have the directories needed to run the alignment pipeline
for this_dir in {logs,error_logs,alignments}; do
    if [ ! -d "${out}"/"${this_dir}" ]
    then
        mkdir -v "${out}"/"${this_dir}"
    else
        echo "Checking whether needed dirs exits: ${out}${this_dir} exists" >> ${output_log}
    fi
done


## -----
## setup and settings
## these usually won't change, but double check:
## be sure the SNP set is two columns for the reference genome
## and that it has the correct line breaks (safe in xcode once, if necessary)
SNPs='/home/albertf/shared/SNPSets/SNPs_Maggie_170809_BY_positions.txt'
head -n 4 ${SNPs}
genome='/home/albertf/shared/genomes/sacCer3.fa'
head -n 4 ${genome}
n_threads=24
PE='PE'
pair_ID1='_R1'
pair_ID2='_R2'
postfix='_001'
file_type='.fastq.gz'


## need to be careful that the underscores in front of ${pair_ID1} match
file_roots=( $(ls ${raw} | grep "${pair_ID1}" | sed "s/${pair_ID1}.*//") )
## make sure you've found your files, e.g.,: 
file_root=${file_roots[10]} #this lists the 10th item in this object
echo ${file_root}


## -----
## find raw files and align
for current_root in ${file_roots[@]}; do #[@] lists all the variables in that object

    file_r1=${raw}${current_root}${pair_ID1}${postfix}${file_type}
    file_r2=${raw}${current_root}${pair_ID2}${postfix}${file_type}

    ## 1. first loop -> ensure we have the first read file
    if [ ! -e ${file_r1} ]
    then
        echo "ERROR: File ${file_r1} not found. Exiting"; exit 1
    else
        echo -e "Found ${file_r1}; now will look for r2\n" >> ${output_log}
    fi

    ## 2. 2nd loop -> ensure we have the second read file
    if [ ! -e ${file_r2} ]
    then
        echo "ERROR: File ${file_r2} not found. Exiting"; exit 1
    else
        echo -e "Found ${file_r2}; starting alignment" >> ${output_log}
    fi

    bwa mem -t ${n_threads} ${genome} ${file_r1} ${file_r2} \
	| samtools sort -@${n_threads} -O BAM -o ${BWA_data}"/"${current_root}_sort.bam -

    ## the above will either produce a sam or exit with an error message
    samtools view -q 30 ${BWA_data}"/"${current_root}_sort.bam \
        | grep 'XS:i:0' \
        | samtools view -b -T ${genome} - \
                   > ${BWA_data}"/"${current_root}_sort_filtered.bam

    ## get rid of PCR duplicates
samtools rmdup -S ${BWA_data}"/"${current_root}_sort_filtered.bam ${BWA_data}"/"${current_root}_sort_filtered_rmdup.bam
# rmdup is obsolete now. man page says us markdup instead

    ## counting the coverage per SNPs
    samtools mpileup -vu -t INFO/AD -l ${SNPs} -f ${genome} ${BWA_data}"/"${current_root}_sort_filtered_rmdup.bam > ${BWA_data}"/"${current_root}_sort_filtered_rmdup.vcf

    echo "done"
done
