#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# optional args
option_list <- list(
   make_option(c("-p", "--prefix"), default = "infile",
               help = "Prefix to use for naming plots (default = derived from infile)"),
   make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
               help = "Verbose mode (default = %default)")
)

# create your own usage
opt_parse <- OptionParser(usage = "%prog [options] dwgsim_eval.output",
                          option_list = option_list)

# set positional_arguments to TRUE
opt <- parse_args(opt_parse, positional_arguments = TRUE)

# print usage if no positional args provided
if (length(opt$args) == 0){
   print_help(opt_parse)
   stop("Please provide an input file")
}

infile <- opt$args[1]

prefix <- vector()
if (grepl("infile", opt$options$prefix)){
   prefix <- strsplit(infile, "\\.")[[1]][1]
} else {
   prefix <- paste0(dirname(infile), '/', opt$options$prefix)
}

suppressPackageStartupMessages(library(tidyverse))

my_cols <- c(
  'mq', 'mc', 'mi', 'mu', 'um', 'uu', 'mt',
  'mct', 'mit', 'mut', 'umt', 'uut', 'mtt',
  'recall', 'ppv', 'fdr',
  'recallt', 'ppvt', 'fdrt'
)

my_col_types <- cols(
  'i', 'i', 'i', 'i', 'i', 'i', 'i',
  'i', 'i', 'i', 'i', 'i', 'i',
  'd', 'd', 'd',
  'd', 'd', 'd'
)

metrics <- read_table(infile, comment = "#", col_names = my_cols, col_types = my_col_types)

if (opt$options$verbose){
   metrics %>%
      summarise(total_read = sum(mt), accuracy = sum(mc)/sum(mt))
}

ggsave(
  paste0(prefix, "_accuracy.png"),
  ggplot(metrics, aes(mq, mc/mt*100)) +
     geom_point() +
     geom_line() +
     ylim(c(0, 100))
)

ggsave(
  paste0(prefix, "_sensitivity.png"),
  ggplot(metrics, aes(mq, recallt)) +
     geom_point() +
     geom_line() +
     ylim(c(0.5, 1))
)

ggsave(
  paste0(prefix, "_ppv.png"),
  ggplot(metrics, aes(mq, ppvt)) +
     geom_point() +
     geom_line()
)

ggsave(
  paste0(prefix, "_fdr.png"),
  ggplot(metrics, aes(mq, fdrt)) +
     geom_point() +
     geom_line()
)

message("Done")
quit()

