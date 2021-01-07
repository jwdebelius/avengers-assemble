#!/usr/bin/env Rscript

# Sets up the libraries we need
suppressMessages(library("argparser"))
suppressMessages(library("vegan"))


# Sets up an argparser to make this tractable
parser <- arg_parser(
  description="Basic Adonis Functionality for age and set"
)
parser <- add_argument(parser,
  '--distance-matrix', 
  help='the distance matrix',
)
parser <- add_argument(parser,
  '--metadata', 
  help=('Path to the metadata file'),
)
parser <- add_argument(parser,
  '--output',
  help=('The location for the summary file')
)

args <- parse_args(parser)

# reads in the distance matrix as a dataframe and sets the rownames as the first column
df <- read.csv2(args$distance_matrix, sep='\t', row.names='X')
map <- read.csv2(args$metadata, sep='\t', row.names='sample.id')
map <- map[row.names(df),]

dm  <- dist(df[row.names(df), row.names(df)])

age.child <- adonis2(dm ~ age, data=map)

res <- data.frame(
  age = c(age$R2[1], age$`Pr(>F)`[1], age$F[1]),
  row.names = c('R2', 'p_999', 'F')
  )

write.table(res, file=args$output, quote=FALSE, sep='\t', col.names = NA)
