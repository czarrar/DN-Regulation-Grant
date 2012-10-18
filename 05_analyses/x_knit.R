#!/home/data/PublicProgram/R/bin/Rscript --vanilla

args <- commandArgs(TRUE)

if (length(args) != 1) {
    cat("x_knit.R script.Rmd\n")
    quit(status=1)
}
    
library(tools)
library(knitr)

opts_knit$set(out.format = "html")
thm <- knit_theme$get('acid')
knit_theme$set(thm)

infile <- file_path_as_absolute(args[1])
#stylesheet <- file_path_as_absolute('foundation/stylesheets/foundation.css')

setwd("/home/data/Projects/CCD/reports")
outfile <- paste("./", file_path_sans_ext(basename(infile)), ".html", sep="")

#knit2html(infile, output=outfile, stylesheet=stylesheet)
knit2html(infile, output=outfile)
