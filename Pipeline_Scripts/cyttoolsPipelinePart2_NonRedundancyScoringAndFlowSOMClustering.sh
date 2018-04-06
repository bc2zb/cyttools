#!/bin/bash

# Fill in the varibles below with the locations of your various files

CYTTOOLS_LOCATION="/Path/To/cyttools/"


WORKING_DIR="/Path/To/Working/Directory/"
FCS_DIR="/Path/To/FCS/Files/"

RESULTS_BLANKS="$WORKING_DIR/cyttoolsBlanksResults/"
RESULTS_CLUSTERING="$WORKING_DIR/cyttoolsClusteringResults/"
RESULTS_DIFFERENTIAL="$WORKING_DIR/cyttoolsDifferentrialResults/"

PANEL="$RESULTS_BLANKS/panelFile.txt"
METADATA="$RESULTS_BLANKS/MetaDataFile.txt"

CLUSTERING_DIR="$RESULTS_CLUSTERING/CLUSTERED_FCS/"
PHENO_DIR="$RESULTS_CLUSTERING/PHENOTYPED_FCS/"
PHENOCODES="$RESULTS_CLUSTERING/PhenoCodes.txt"

cd $CYTTOOLS_LOCATION

Rscript cyttools.R --computeNRS "$FCS_DIR" "$PANEL" "$RESULTS_BLANKS"

NRS_PANEL="$RESULTS_BLANKS/nrsPanelFile.txt"

Rscript cyttools.R --cluster=FlowSOM "$FCS_DIR" "$NRS_PANEL" "$RESULTS_CLUSTERING"
Rscript cyttools.R --cluster=FlowType "$FCS_DIR" "$NRS_PANEL" "$RESULTS_CLUSTERING"
Rscript cyttools.R --phenoConsensusClustering "$CLUSTERING_DIR" "$PHENO_DIR" "$NRS_PANEL" "$RESULTS_CLUSTERING"
