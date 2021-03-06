#!/bin/bash

# Fill in the varibles below with the locations of your various files

CYTTOOLS_LOCATION="/Path/To/cyttools/"

WORKING_DIR="/Path/To/Working/Directory/"
FCS_DIR="/Path/To/FCS/Files/"

RESULTS_BLANKS="$WORKING_DIR/cyttoolsBlanksResults/"
RESULTS_CLUSTERING="$WORKING_DIR/cyttoolsClusteringResults/"
RESULTS_DIFFERENTIAL="$WORKING_DIR/cyttoolsDifferentrialResults/"

PANEL="$RESULTS_BLANKS/nrsPanelFile.txt"
METADATA="$RESULTS_BLANKS/MetaDataFile.txt"

BATCH_RESULTS_DIR="$RESULTS_CLUSTERING/BatchFlowTypeDataPrepWorkspaces/"
BATCH_DIR="$BATCH_RESULTS_DIR/BatchFlowTypeResults/"

cd $CYTTOOLS_LOCATION

Rscript cyttools.R --cluster=BatchFlowTypeDataPrep "$FCS_DIR" "$PANEL" "$RESULTS_CLUSTERING"

cd $BATCH_RESULTS_DIR

for i in $(find $BATCH_RESULTS_DIR -type f -d 1); do
  cd $CYTTOOLS_LOCATION

  Rscript cyttools.R --batchFlowType "$i" "$BATCH_DIR"

  echo $i
done

# perfrom clustering analysis, WARNING FlowType takes a long time to run and will eat up most of your memory

Rscript cyttools.R --batchFlowTypeDataMerge "$FCS_DIR" "$PANEL" "$RESULTS_CLUSTERING" "$BATCH_DIR"

ABNDNC_FILE_NAME="nodeAbundanceFeatureTable.txt"
ABNDNC_FEATURE_TABLE="$RESULTS_CLUSTERING/$ABNDNC_FILE_NAME"

Rscript cyttools.R --compDiffAbndnc  "$PANEL" "$ABNDNC_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"

EXPR_FILE_NAME="nodeExpressionFeatureTable.txt"
EXPR_FEATURE_TABLE="$RESULTS_CLUSTERING/$EXPR_FILE_NAME"

Rscript cyttools.R --compDiffExpr "$PANEL" "$EXPR_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"

                                                                                                                         