library(plyr)
library(ggplot2)
library(RColorBrewer)
library(robustbase)
basedir <- "/home2/data/Projects/CCD"
scriptdir <- file.path(basedir, "scripts")
datadir <- file.path(basedir, "scripts/data")

# Get design files
## CCB
ccb_coherent   <- read.table(file.path(scriptdir, "04_msit_task/CCB_coherent.1D"))[,1] == 1
ccb_coherent <- ccb_coherent[1:(length(ccb_coherent)-2)]
ccb_incoherent <- read.table(file.path(scriptdir, "04_msit_task/CCB_incoherent.1D"))[,1] == 1
ccb_incoherent <- ccb_incoherent[1:(length(ccb_incoherent)-2)]
## CCD
ccd_coherent   <- read.table(file.path(scriptdir, "04_msit_task/CCD_coherent.1D"))[,1] == 1
ccd_coherent <- ccd_coherent[1:(length(ccd_coherent)-2)]
ccd_incoherent <- read.table(file.path(scriptdir, "04_msit_task/CCD_incoherent.1D"))[,1] == 1
ccd_incoherent <- ccd_incoherent[1:(length(ccd_incoherent)-2)]

# Get subject info
subinfo <- read.csv(file.path(scriptdir, "04_msit_task/z_predesign.csv"))
subinfo$study <- factor(subinfo$study, labels=c("CCB", "CCD"))

# Get the filenames
fnames_df <- ddply(subinfo, .(subject, study, scan, run), function(sdf) {
    if (sdf$study == "CCB") {
        if (sdf$scan == 1) {
            find_path <- file.path(
                            basedir, "CCB/extras/logfiles", 
                            paste(sdf$subject, "_msit_", sdf$run, "_behavior.1D", sep="")
                        )
        } else {
            find_path <- file.path(
                            basedir, "CCB/extras/logfiles", 
                            paste(sdf$subject, "_Scan", sdf$scan, "_msit_", sdf$run, "_behavior.1D", sep="")
                        )
        }
        
    } else if (sdf$study == "CCD") {
        find_path <- file.path(
                        basedir, "behavior/logfiles", 
                        paste(sdf$subject, "_log_MSIT", sdf$scan, "_behav.txt", sep="")
                    )
    } else {
        print(sdf)
        stop("Unknown study ", sdf$study)
    }
    paths <- Sys.glob(find_path)
    if (length(paths) != 1) {
        print(find_path)
        print(paths)
        stop("Error 0 or 2+ paths found")
    }
    c(infiles=paths)
})

# Get time-serieses
tss <- ddply(fnames_df, .(subject, study, scan, run), function(sdf) {
    behav <- as.matrix(read.table(sdf$infiles))
    bad <- behav[,2] == -1 | behav[,3] == 1
    if (sdf$study == "CCB") {
        cc <- behav[ccb_coherent & !bad, 2]
        ii <- behav[ccb_incoherent & !bad, 2]
    } else {
        cc <- behav[ccd_coherent & !bad, 2]
        ii <- behav[ccd_incoherent & !bad, 2]            
    }
    nc <- length(cc)
    ni <- length(ii)
    
    data.frame(
        condition = rep(c("Coherent", "Incoherent"), times=c(nc, ni)), 
        rt = c(cc, ii)
    )
})

