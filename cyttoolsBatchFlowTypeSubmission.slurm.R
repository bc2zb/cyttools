#!/bin/bash
#SBATCH --ntasks=1
#SBATCH -t 24:00:00
#SBATCH -p economy
#SBATCH -A computational_cyto
#SBATCH --mem-per-cpu=16000
#SBATCH --output=result_%a.out

# load R
module load R/openmpi/3.1.1

# Fill in the varibles below with the locations of your various files

CYTTOOLS_LOCATION="/home/bc2zb/cyttools-0.2.1/"
FCS_DIR="/scratch/bc2zb/JXSTUFF/ALL_FILES/AllFiles/"
PANEL="/scratch/bc2zb/JXSTUFF/panelForCyttoolsRun.txt"
METADATA="/scratch/bc2zb/JXSTUFF/MetaDataFile.txt"
RESULTS_BLANKS="/scratch/bc2zb/JXSTUFF/ALL_FILES/cyttoolsBlanksResults/"
RESULTS_CLUSTERING="/scratch/bc2zb/JXSTUFF/ALL_FILES/cyttoolsClusteringResults/"
RESULTS_DIFFERENTIAL="/scratch/bc2zb/JXSTUFF/ALL_FILES/cyttoolsDifferentialResults/"

cd $CYTTOOLS_LOCATION

Rscript cyttools.R --cluster=BatchFlowType "$FCS_DIR" "$PANEL" "$RESULTS_CLUSTERING"


