#!/usr/bin/env Rscript 

library(docopt)
"Usage:
my_program.R [options]

Options:
--Threshold=<Threshold>       [default: 250]
--output=OUTPUT               [default: out.txt]
--color=<color>               [default: FALSE]
" -> doc

opt <- docopt(doc)
print(opts) 