## cyttools

cyttools is an open source pipeline framework for the analysis of high parameter cytometry data. It was presented as a poster at ABRF 2017 and CYTO 2017 and won outstanding posters awards at both conferences. It is currently ready to be used but is undergoing active development.

## Resources for cytometry data analysis

[Introduction to R](http://briancapaldo.com/SlideDeck.html)

[Comparison of Feature Selection Methods](http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23030/full)

cyttools approach to data anlysis is heavily influenced by the paper below.

[F1000 paper on high parameter cytometry analysis](https://f1000research.com/articles/6-748/v1)

## Installing cyttools

Clone the repository to your local machine. Install R from [CRAN](https://cran.r-project.org/)

Next, run the "Install required packages_cyttools.R" script using R. If successful, cyttools should be fully functional.

cyttoolsPipeline.sh contains a generic analysis pipeline. Replace the file and directory paths at the top of the file with those on your system, and you should be able to run the analysis by running the script from the command line.

## Methods

Mass cytometry data is inverse hyberbolic sin transformed using a cofactor of 0.25. [FlowSOM](http://bioconductor.org/packages/release/bioc/html/FlowSOM.html) is used to construct a self-organizing map with a number of grid points equal to the square of the number of lineage markers. Each cell is assigned a phenocode for every lineage marker using [flowType](http://bioconductor.org/packages/release/bioc/html/flowType.html). Each grid point is then immunophenotyped using the aggregated phenocodes of the cells assigned to the grid point. For any given grid point to be assigned an immunophenotype for a particular marker (i.e. CD45+), 75 percent of the cells assigned to the gridpoint must be labelled with the same phenocode for the particular marker. For each immunophenotype observed, number of cells are tabulated to form a hierarchical count table. Every level of the hierarchy is tested for differential abundance between conditions using edgeR with a quasi-likelihood framework as specified by the [cydar](http://bioconductor.org/packages/release/bioc/html/cydar.html) package. 

## References

Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2017). FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.
  http://www.r-project.org, http://dambi.ugent.be.
  
Nima Aghaeepour, Kieran O'Neill and Adrin Jalali (2014). flowType: Phenotyping Flow Cytometry Assays. R package version 2.16.0.

Aaron Lun (2017). cydar: Using Mass Cytometry for Differential Abundance Analyses. R package version 1.2.1.