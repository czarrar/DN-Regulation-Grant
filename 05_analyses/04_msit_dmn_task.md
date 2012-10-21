# Association between MSIT task design and DN signal

Note that actual code is loaded from a different file.


```r
read_chunk("04_msit_dmn_task.R")
```


## Setup


```r
library(plyr)
library(e1071)
library(ggplot2)
library(RColorBrewer)
library(robustbase)
library(reshape)
basedir <- dirname(dirname(getwd()))  # assume running in current direcotry
scriptdir <- file.path(basedir, "scripts/04_msit_task")
datadir <- file.path(basedir, "scripts/data")
oldtheme <- theme_set(theme_bw())
```



```r
network_names <- c("medial visual", "occipital pole visual", "lateral visual", 
    "default network", "cerebellum", "sensorimotor", "auditory", "executive control", 
    "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- 4
tps <- 8:10
```



```r
# only CCD subjects with MSIT
phenos <- read.csv(file.path(datadir, "ccd_totals_touse.csv"))[1:9, -1]
subinfo <- read.csv(file.path(scriptdir, "z_predesign.csv"))
subinfo$study <- factor(subinfo$study, labels = c("CCB", "CCD"))
subinfo$scan <- factor(subinfo$scan)
subinfo$run <- factor(subinfo$run)
subinfo$sex <- factor(subinfo$sex)
# Plot Age
ggplot(subinfo, aes(x = age, fill = ..count..)) + geom_histogram(binwidth = 5) + 
    facet_grid(study ~ .)
```

![plot of chunk subject-info](figure/subject-info1.png) 

```r
# Plot Sex
ggplot(subinfo, aes(x = sex, fill = study)) + geom_bar()
```

![plot of chunk subject-info](figure/subject-info2.png) 


Note that negative signal here are the congruent trials and positive signal are the incongruent trials. What do you think would happen if the positives and negatives got together?


```r
# CCB
ccb_waver <- read.table(file.path(scriptdir, "level1_ccb_template.mat"), skip = 5)
ccb_waver <- as.matrix(ccb_waver)[, c(1, 3)]
ccb_waver <- ccb_waver[, 2] - ccb_waver[, 1]
tmpdf <- data.frame(time = seq(0, by = 1.75, length.out = length(ccb_waver)), 
    predicted_signal = ccb_waver)
ggplot(tmpdf, aes(time, predicted_signal)) + geom_line() + xlab("Time (secs)") + 
    ylab("Predicted Signal") + ggtitle("For CCB Subjects")
```

![plot of chunk task-waver](figure/task-waver1.png) 

```r
# CCD
ccd_waver <- read.table(file.path(scriptdir, "level1_ccd_template.mat"), skip = 5)
ccd_waver <- as.matrix(ccd_waver)[, c(1, 3)]
ccd_waver <- ccd_waver[, 2] - ccd_waver[, 1]
tmpdf <- data.frame(time = seq(0, by = 2, length.out = length(ccd_waver)), predicted_signal = ccd_waver)
ggplot(tmpdf, aes(time, predicted_signal)) + geom_line() + xlab("Time (secs)") + 
    ylab("Predicted Signal") + ggtitle("For CCD Subjects")
```

![plot of chunk task-waver](figure/task-waver2.png) 



```r
load(file.path(datadir, "ccb+ccd_rts_all.rda"))  # rts
rts$scan <- factor(rts$scan)
rts$run <- factor(rts$run)
```



```r
# this loads the 'tss' object with attr(tss, 'split_labels') to get how
# stuff should be organized
load(file.path(basedir, "scripts/data/ccb+ccd_time_series.rda"))
splitter <- attr(tss, "split_labels")
splitter$index <- 1:nrow(splitter)
# only look at time-series for MSIT with associated RT info note: there's
# proly a better
splitter <- ddply(splitter, .(study, subject, condition, scan, run), function(x) {
    if (x$condition == "REST") 
        return(x)
    has_any <- any(ddply(subset(rts, trial == 1), .(study, subject, condition, 
        scan, run), function(y) {
        as.character(x$subject) == as.character(y$subject) & x$scan == y$scan & 
            x$run == y$run
    })$V1)
    if (has_any) 
        return(x) else return(data.frame())
})
```


## Group Average

Get the xmin and xmax for showing blocks.


```r
# CCB
ccb_coherent <- read.table(file.path(scriptdir, "CCB_coherent.1D"))[, 1]
ccb_incoherent <- read.table(file.path(scriptdir, "CCB_incoherent.1D"))[, 1]
tmp <- diff(ccb_coherent - ccb_incoherent)
# +1 to adjust for relative; +3 to adjust for HRF delay
ccb_event_tpts <- data.frame(xmin = c(0, which(tmp != 0) + 1 + 3), xmax = c(which(tmp != 
    0) + 1 + 3, length(ccb_coherent)), block = factor(c("Fixation", rep(c("Coherent", 
    "Incoherent"), length.out = sum(tmp != 0) - 1), "Fixation")))
ccb_event_tpts$xmin <- ccb_event_tpts$xmin * 1.75
ccb_event_tpts$xmax <- ccb_event_tpts$xmax * 1.75
# CCD
ccd_coherent <- read.table(file.path(scriptdir, "CCD_coherent.1D"))[, 1]
ccd_incoherent <- read.table(file.path(scriptdir, "CCD_incoherent.1D"))[, 1]
tmp <- diff(ccd_coherent - ccd_incoherent)
# +1 to adjust for relative; +3 to adjust for HRF delay
ccd_event_tpts <- data.frame(xmin = c(0, which(tmp != 0) + 1 + 3), xmax = c(which(tmp != 
    0) + 1 + 3, length(ccd_coherent)), block = factor(c("Fixation", rep(c("Coherent", 
    "Incoherent"), length.out = sum(tmp != 0) - 1), "Fixation")))
ccd_event_tpts$xmin <- ccd_event_tpts$xmin * 2
ccd_event_tpts$xmax <- ccd_event_tpts$xmax * 2
```


### BOLD Time Series

Compute the group average MSIT BOLD signal.









































