#!/bin/bash

#### WARNING ####

# The most common parameters used for setting up scripts for running jobs on SLURM queueing system are reported
# The SLURM scripts have to be intended as bash script, hence the interpreter always occupies the first line.
# Single hash (#) is used either for commenting and for introducing variable (#SBATCH), in the form of a flag, that will
# be interpreted by the sbatch manager, i.e. running and queuing the script "scriptname.sh" as
#
# sbatch scriptname.sh
#
# The double hash is used for commenting the variables in order to skip their usage (e.g. in the present case this
# can be referred to the use of a parallel environment or the multiple task issue).

#### PARTITION WHERE DIRECT YOUR JOB/S
#SBATCH --partition=batch

#### GENERIC PARAMETERS ####

# HARD TIME LIMIT
# Specify hard time limit for the job (currently, there is no hard limit imposed by administrators on this machine)
#SBATCH --time=24:00:00

# MAIL ALERT ## Not working right now
# Send an email when the name.surname@ifom.eujob finishes or if it is aborted (by default no email is sent).
# The possible values are – when the job begins (begin), ends (end), is fail (fail) or never (none) – default.
##SBATCH --mail-type=begin
##SBATCH --mail-user=edoardo.bellini@ifom.eu

# JOB NAME
# Give job a name (possibly avoid spaces)
#SBATCH --job-name=rna_seq

# ERROR FILE
# Error file name (possibly avoid spaces)
#SBATCH --error="error.err"


# OUTPUT FILE
# Output file (possibly avoid spaces)
#SBATCH --output="output.out"

# UNIQUE OUTPUT FILE
# Combine output and error files into a single file (indicated using the -o flag) could be achived just not specifing "--error" (see few line above)

# WORKING DIRECTORY
# By default a SLURM batch script is exectued in the directory where it has been launched

# CORES AND MEMORY
# Resources required for the submitted job. Request a strict series of resources only if really needed, otherwise the queuing system could hold your job submission for a very long and unwanted time (e.g., be careful with the requested memory, here commented).
##SBATCH --cpus-per-task=1 
##SBATCH --mem-per-cpu=4G  
#SBATCH --mem=300GB
#SBATCH --ntasks=50

# JOB ARRAYS
# Submit an array job with more than one task (e.g., 4 tasks)
##SBATCH --array=0-4
#
# Then use the $SLURM_ARRAY_TASK_ID environment variable to select the appropriate input file from bash array
# Bash array index starts from 0, so we need to subtract one from $SLURM_ARRAY_TASK_ID value
#inputs=(input1 input2 input3)
#index=$(($SLURM_ARRAY_TASK_ID-1))
#taskinput=${inputs[$index]}
#
# and use it like this
#./my_program $taskinput   result.$SLURM_ARRAY_TASK_ID

# INFORMATION AND LOG FILE
# Keep track of information related to the current job in the log file (see the previous point)
# The batch environment variables are used

echo "=========================================================="
echo Running in the directory `pwd`
echo Running on node $hostname
echo "Start date : $(date)"
echo "Job name : $SBATCH_JOB_NAME"
echo "Job ID : $SLURM_JOB_ID  $SLURM_ARRAY_TASK_ID"
echo "Number of slots requested : $SLURM_NPROCS"
echo "=========================================================="

# COMMAND
# Program name or command and its options and arguments
# The last line of the script MUST be left empty, otherwise the system will detect an error in processing it
#
# This is just for testing you are really doing something in the place the process is supposed to run,
# i.e. a new empty file is written
touch empty.log

# This is the test command, i.e. a straight binary executable the reads two input file, uses 4 cores and generates many files
# The run takes normally less than a minute


nextflow run nf-core/rnaseq \
	--outdir output \
	-profile slurm,singularity \
	-c /storage/data/l-hpc/SLURM_tests/nf-core/CONFIG/conf/ifom.config \
	--input /storage/data/GS/data/genewiz_2024-08-06/my_samplesheet_test.csv \
	--fasta /storage/data/GS/data/genome/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--gtf /storage/data/GS/data/genome/index/Homo_sapiens.GRCh38.104.gtf \
	-r 3.13.2 \
	--star_index /storage/data/GS/data/genome/index/star/ \
	--salmon_index /storage/data/GS/data/genome/index/salmon \
	--rsem_index /storage/data/GS/data/genome/rsem