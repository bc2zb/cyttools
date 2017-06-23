#!/usr/bin/env Rscript

require(docopt)
require(methods)

"
Usage:
rand-unif [--number=<number>] [--min=<min>] [--max=<max>]
rand-unif (-h | --help | --version)

Description:   This program generates random unifom numbers.
Options:
--version       Show the current version.
--number=<num>  [default: 1] The number of random numbers to generate.
--min=<min>     [default: 0] The lowest value a random number can have.
--max=<max>     [default: 1] The highest value a random number can have.
" -> doc


args    <- docopt(doc)

minimum <- as.numeric(args $ `--min`)
maximum <- as.numeric(args $ `--max`)
number  <- as.numeric(args $ `--number`)

cat(paste0(runif(number, minimum, maximum), collapse = '\n'), '\n')