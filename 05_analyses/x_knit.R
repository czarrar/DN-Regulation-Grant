#!/usr/bin/env Rscript --vanilla

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

d <- dirname(args[1])
infile <- file_path_as_absolute(args[1])
if (file_ext(infile) != "Rmd") stop("Input file must have extension .Rmd")
#stylesheet <- file_path_as_absolute('foundation/stylesheets/foundation.css')

argv <- commandArgs(trailingOnly = FALSE)
reportdir <- file.path(dirname(dirname(dirname(file_path_as_absolute(substring(argv[grep("--file=", argv)], 8))))), "reports")
print(reportdir)
setwd(reportdir)
if (d != "." && !file.exists(d))
    dir.create(d, recursive=TRUE)
outfile <- paste(d, "/", file_path_sans_ext(basename(infile)), ".html", sep="")

#knit2html(infile, output=outfile, stylesheet=stylesheet)
knit2html(infile, output=outfile)
