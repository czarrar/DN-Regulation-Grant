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
## Get condition labels
### CCB
ccb_trials                                  = ccb_coherent | ccb_incoherent
ccb_conditions                              = vector("character", sum(ccb_trials))
ccb_conditions[ccb_coherent[ccb_trials]]    = "Coherent"
ccb_conditions[ccb_incoherent[ccb_trials]]  = "Incoherent"
### CCD
ccd_trials                                  = ccd_coherent | ccd_incoherent
ccd_conditions                              = vector("character", sum(ccd_trials))
ccd_conditions[ccd_coherent[ccd_trials]]    = "Coherent"
ccd_conditions[ccd_incoherent[ccd_trials]]  = "Incoherent"

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
                        paste(sdf$subject, "_log_MSIT", sdf$run, "_behav.txt", sep="")
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

# Get time-serieses (don't include trials with errors and trials with RT < 200)
rts <- ddply(fnames_df, .(subject, study, scan, run), function(sdf) {
    if (sdf$study == "CCB") {
        behav <- as.matrix(read.table(sdf$infiles))
        bad <- behav[,2] == -1 | behav[,3] == 1 | behav[,2] < 200
        cc <- behav[ccb_coherent & !bad, 2]
        ii <- behav[ccb_incoherent & !bad, 2]
    } else {
        behav <- as.matrix(read.table(sdf$infiles, skip=2))
        bad <- behav[,2] == -1 | behav[,3] == 1 | behav[,2] < 200
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
# Exclude subjects with no responses or errors > 10
nccb <- sum(ccb_coherent+ccb_incoherent)
nccd <- sum(ccd_coherent+ccd_incoherent)
rts <- ddply(rts, .(study,subject,scan,run), function(x) {
    if (x$study[1]=="CCB") nbad <- nccb-nrow(x)
    else if (x$study[1]=="CCD") nbad <- nccd-nrow(x)
    if (nbad > 10) return(data.frame())
    else return(x)
})    
rts$subject <- factor(rts$subject)
# Save
save(rts, file=file.path(datadir, "ccb+ccd_rts.rda"))

# Get time-serieses (include trials with errors as NA and trials with RT < 200)
rts <- ddply(fnames_df, .(subject, study, scan, run), function(sdf) {
    if (sdf$study == "CCB") {
        behav <- as.matrix(read.table(sdf$infiles))
        bad <- behav[,2] == -1 | behav[,3] == 1 | behav[,2] < 200
        behav[bad,2] <- NA
        rt <- behav[ccb_coherent|ccb_incoherent, 2]
        condition <- ccb_conditions
    } else {
        behav <- as.matrix(read.table(sdf$infiles, skip=2))
        bad <- behav[,2] == -1 | behav[,3] == 1 | behav[,2] < 200
        behav[bad,2] <- NA
        rt <- behav[ccd_coherent|ccd_incoherent,2]
        condition <- ccd_conditions
    }
    
    data.frame(
        trial = 1:length(condition), 
        condition = condition, 
        rt = rt
    )
})
# Exclude subjects with no responses or errors > 10
nccb <- sum(ccb_coherent+ccb_incoherent)
nccd <- sum(ccd_coherent+ccd_incoherent)
rts <- ddply(rts, .(study,subject,scan,run), function(x) {
    nbad <- sum(is.na(x$rt))
    if (nbad > 10) return(data.frame())
    else return(x)
})
rts$subject <- factor(rts$subject)
# Save
save(rts, file=file.path(datadir, "ccb+ccd_rts_all.rda"))
