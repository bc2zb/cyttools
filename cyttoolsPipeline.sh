#!/bin/bash

# Fill in the varibles below with the locations of your various files

CYTTOOLS_LOCATION="/PATH/TO/CYTTOOLS/DIRECTORY/"
FCS_DIR="/PATH/TO/FCS/FILES/"
PANEL="/PATH/TO/PANEL/FILE.TXT"
METADATA="/PATH/TO/METADATA/FILE.TXT"
RESULTS_BLANKS="/PATH/TO/CLUSTERING/BLANKS/DIRECTORY/"
RESULTS_CLUSTERING="/PATH/TO/CLUSTERING/RESULTS/DIRECTORY/"
RESULTS_DIFFERENTIAL="/PATH/TO/DIFFERENTIAL/RESULTS/DIRECTORY/"

cd $CYTTOOLS_LOCATION

# if needed, generate panel and meta data files. DO THIS BEFORE YOU RUN THIS SCRIPT!

# Rscript cyttools.R --makePanelBlank "$FCS_DIR" "$RESULTS_BLANKS"
# Rscript cyttools.R --makeMetaDataBlank "$FCS_DIR" "$RESULTS_BLANKS"

# perfrom clustering analysis, WARNING FlowType takes a long time to run and will eat up most of your memory

Rscript cyttools.R --cluster=FlowSOM "$FCS_DIR" "$PANEL" "$RESULTS_CLUSTERING"

# perform differential analysis, Rscript cyttools.R --cluster will automatically generate feature tables to be used in this command, and the above BASH command moves them to where the next command can find them

ABNDNC_FILE_NAME="nodeAbundanceFeatureTable.txt"
ABNDNC_FEATURE_TABLE="$RESULTS_CLUSTERING/$ABNDNC_FILE_NAME"

Rscript cyttools.R --compDiffAbndnc  "$PANEL" "$ABNDNC_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"

EXPR_FILE_NAME="nodeExpressionFeatureTable.txt"
EXPR_FEATURE_TABLE="$RESULTS_CLUSTERING/$EXPR_FILE_NAME"

Rscript cyttools.R --compDiffExpr "$PANEL" "$EXPR_FEATURE_TABLE" "$METADATA" "$RESULTS_DIFFERENTIAL"
