## cyttools

cyttools is an open source pipeline framework for the analysis of high parameter cytometry data. It was presented as a poster at ABRF 2017 and CYTO 2017 and won outstanding posters awards at both conferences. It is currently ready to be used but is undergoing active development.

## Resources for cytometry data analysis

[Introduction to R](http://briancapaldo.com/SlideDeck.html)

[Comparison of Feature Selection Methods](http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23030/full)

[Author's GitHub](https://github.com/lmweber)

cyttools approach to data anlysis is heavily influenced by the paper below.

[F1000 paper on high parameter cytometry analysis](https://f1000research.com/articles/6-748/v1)

## Installing cyttools

Clone the repository to your local machine. Install R from [CRAN](https://cran.r-project.org/)

Next, run the "Install required packages_cyttools.R" script using R. If successful, cyttools should be fully functional.

cyttoolsPipeline.sh contains a generic analysis pipeline. Replace the file and directory paths at the top of the file with those on your system, and you should be able to run the analysis by running the script from the command line.

Mass cytometry data was inverse hyberbolic sin transformed using a cofactor of 0.25. Lineage markers were used to construct a self-organizing map with 529 grid points. Phenocodes for every cell were derived using flowType and each grid point was immunophenotyped using the phenocodes of the cells assigned to the grid point. For each marker, grid points containing at least 75 percent of the same cell type were labeled as those cell type. For each phenotype observed, number of cells were tabulated to form a hierarchical count table. Every level of the hierarchy was tested for differential abundance between conditions using edgeR with a quasi likelihood framework. 