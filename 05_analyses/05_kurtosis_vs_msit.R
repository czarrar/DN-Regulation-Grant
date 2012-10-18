## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
basedir <- "/home2/data/Projects/CCD"

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10

## @knitr -----------break-------------


## @knitr phenotypes
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)[1:9,] # only for 9 subjects with MSIT

# Read in Task Waver
task_waver <- read.table(file.path(basedir, "rois/waver_msit_design.1D"))
task_waver <- as.vector(as.matrix(task_waver))

## @knitr timeseries
# Run 1
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_01/rsn10.1D")))
msit_run1_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
msit_run1_tcs <- msit_run1_tcs[,,dmn]
# Run 2
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_02/rsn10.1D")))
msit_run2_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
msit_run2_tcs <- msit_run2_tcs[,,dmn]
