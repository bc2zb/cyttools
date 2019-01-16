########################################################################################
## Install required packages for Advanced R course section from CRAN and Bioconductor
required_packages_CRAN <- c("tidyverse", "reshape2", "docopt", "methods", "parallel", "devtools")
required_packages_BioC <- c("flowCore", "limma", "FlowSOM", "flowType", "ncdfFlowSet")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(required_packages_CRAN)

BiocManager::install(required_packages_BioC)

########################################################################################
## Tests whether all required packages 
## Press the 'Source' button in the top right corner of this pane and check 
## whether the output in the Console pane confirms that all packages are installed

required_packages <- c(required_packages_CRAN, required_packages_BioC)

installed_packages <- required_packages %in% installed.packages()[,"Package"]
missing_packages <- required_packages[!installed_packages]
if ( length(missing_packages) > 0 ) {
	cat(sprintf('FOLLOWING PACKAGES NEED TO BE INSTALLED STILL:\n\t%s\n',
		paste(missing_packages, collapse=', ')))
} else{
	cat('ALL PACKAGES ARE INSTALLED, cyttools should work')
}

