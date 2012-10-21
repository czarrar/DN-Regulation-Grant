# This script is intended to gather CCB and CCD subject network time-series into
# one rda. It's going to be a merry christmas when this script is done.

library(plyr)
basedir <- "/home2/data/Projects/CCD"
scriptdir <- file.path(basedir, "scripts/05_analyses")
datadir <- file.path(basedir, "scripts/data")

# Read in subject info
subinfo <- read.csv(file.path(datadir, "ccb+ccd_filtered_all.csv"))
subinfo$scan <- factor(subinfo$scan)
subinfo$run <- factor(subinfo$run)

# Get data frame with associated TS filepaths
fnames_df <- ddply(subinfo, .(subject, study, condition, scan, run), function(sdf) {
    if (sdf$study == "CCB") {
        find_path <- file.path(
                        basedir, "CCB/analysis", 
                        paste(sdf$subject, "_scan", sdf$scan, "_", sdf$condition, sdf$run, "_*_RSN_TCs.1D", sep="")
                    )
    } else if (sdf$study == "CCD") {
        find_path <- file.path(
                        basedir, "analysis/subjects", sdf$subject, 
                        tolower(sdf$condition), sprintf("run_%02i", sdf$run), 
                        "rsn10.1D"
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
tss <- dlply(fnames_df, .(subject, study, condition, scan, run), function(sdf) {
    as.matrix(read.table(sdf$infiles))
})
## can use attr(tss, "split_labels") to get how the splitting was done

save(tss, file=file.path(datadir, "ccb+ccd_time_series_all.rda"))
