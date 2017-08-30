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


cd $CYTTOOLS_LOCATION

Rscript cyttools.R --computeNRS "$FCS_DIR" "$PANEL" "$RESULTS_BLANKS"

NRS_PANEL="$RESULTS_BLANKS/nrsPanelFile.txt"

Rscript cyttools.R --cluster=FlowSOM "$FCS_DIR" "$NRS_PANEL" "$RESULTS_CLUSTERING"

ABNDNC_FILE_NAME="nodeAbundanceFeatureTable.txt"
ABNDNC_FEATURE_TABLE="$RESULTS_CLUSTERING/$ABNDNC_FILE_NAME"

Rscript cyttools.R --compDiffAbndnc  "$NRS_PANEL" "$ABNDNC_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"

EXPR_FILE_NAME="nodeExpressionFeatureTable.txt"
EXPR_FEATURE_TABLE="$RESULTS_CLUSTERING/$EXPR_FILE_NAME"

Rscript cyttools.R --compDiffExpr "$NRS_PANEL" "$EXPR_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 