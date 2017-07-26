#!/bin/bash

# Fill in the varibles below with the locations of your various files

CYTTOOLS_LOCATION="/PATH/TO/CYTTOOLS/DIRECTORY/""
FCS_DIR="/PATH/TO/FCS/FILES/"
PANEL="/PATH/TO/PANEL/FILE.TXT"
METADATA="/PATH/TO/METADATA/FILE.TXT"
RESULTS_CLUSTERING="/PATH/TO/CLUSTERING/RESULTS/DIRECTORY/"
RESULTS_DIFFERENTIAL="/PATH/TO/DIFFERENTIAL/RESULTS/DIRECTORY/"

cd $CYTTOOLS_LOCATION

# if needed, generate panel and meta data files. DO THIS BEFORE YOU RUN THIS SCRIPT!

# Rscript cyttools.R --makePanelBlank "$FCS_DIR"
# Rscript cyttools.R --makeMetaDataBlank "$FCS_DIR"

# perfrom clustering analysis, WARNING FlowType takes a long time to run and will eat up most of your memory

Rscript cyttools.R --cluster=FlowSOM "$FCS_DIR" "$PANEL"

# move clustering results directory from cyttoolsResults to your experiment location, this script will automatically create the directory, you should not make it beforehand

mkdir "$RESULTS_CLUSTERING"
find cyttoolsResults/* -type f -print | xargs -I {} mv {} "$RESULTS_CLUSTERING"

# perform differential analysis, Rscript cyttools.R --cluster will automatically generate feature tables to be used in this command, and the above BASH command moves them to where the next command can find them

ABNDNC_FEATURE_TABLE='$RESULTS_CLUSTERINGnodeAbundanceFeatureTable.txt'
EXPR_FEATURE_TABLE='$RESULTS_CLUSTERINGnodeExpressionFeatureTable.txt'

Rscript cyttools.R --compDiffAbndnc "$PANEL" "$ABNDNC_FEATURE_TABLE" "$METADATA"
Rscript cyttools.R --compDiffExpr "$PANEL" "$EXPR_FEATURE_TABLE" "$METADATA"

# move differential analysis results directory from cyttoolsResults to your experiment location

mkdir "$RESULTS_DIFFERENTIAL"
find cyttoolsResults/* -type f -print | xargs -I {} mv {} "$RESULTS_DIFFERENTIAL"

