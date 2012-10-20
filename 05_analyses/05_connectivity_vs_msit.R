## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
library(robustbase)
library(RColorBrewer)
basedir <- dirname(dirname(getwd()))    # assume running in current direcotry
scriptdir <- file.path(basedir, "scripts")
datadir <- file.path(basedir, "scripts/data")
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
orig_network_names <- c("Medial Visual", "Occipital Pole Visual", "Lateral Visual", "Default Network", "Cerebellum", "Sensorimotor", "Auditory", "Executive Control", "Right Frontoparietal", "Left Frontoparietal")
network_names <- gsub(" ", ".", tolower(orig_network_names))
dmn <- which(network_names == "default.network")
tps <- 8:10


## @knitr -----------break-------------


###
# Setup
###

## @knitr subject-info
fname <- file.path(scriptdir, "04_msit_task/z_predesign.csv")
subinfo <- read.csv(fname)
subinfo$study <- factor(subinfo$study, labels=c("CCB", "CCD"))
subinfo$sex <- factor(subinfo$sex)
# Plot Age
ggplot(subinfo, aes(x=age, fill=..count..)) + geom_histogram(binwidth=5) + facet_grid(study ~ .)
# Plot Sex
ggplot(subinfo, aes(x=sex, fill=study)) + geom_bar()

## @knitr task-waver
# CCB
ccb_waver <- read.table(file.path(scriptdir, "04_msit_task/level1_ccb_template.mat"), skip=5)
ccb_waver <- as.matrix(ccb_waver)[,c(1,3)]
ccb_waver <- ccb_waver[,2] - ccb_waver[,1]
tmpdf <- data.frame(
    time = seq(0, by=1.75, length.out=length(ccb_waver)), 
    predicted_signal = ccb_waver
)
ggplot(tmpdf, aes(time, predicted_signal)) + 
    geom_line() + 
    xlab("Time (secs)") + 
    ylab("Predicted Signal") +
    ggtitle("For CCB Subjects")
# CCD
ccd_waver <- read.table(file.path(scriptdir, "04_msit_task/level1_ccd_template.mat"), skip=5)
ccd_waver <- as.matrix(ccd_waver)[,c(1,3)]
ccd_waver <- ccd_waver[,2] - ccd_waver[,1]
tmpdf <- data.frame(
    time = seq(0, by=2, length.out=length(ccd_waver)), 
    predicted_signal = ccd_waver
)
ggplot(tmpdf, aes(time, predicted_signal)) + 
    geom_line() + 
    xlab("Time (secs)") + 
    ylab("Predicted Signal") +
    ggtitle("For CCD Subjects")

## @knitr msit-ts
# this loads the 'tss' object with attr(tss, 'split_labels') to get how stuff should be organized
load(file.path(datadir, "ccb+ccd_time_series.rda"))
splitter <- attr(tss, 'split_labels')
splitter$index <- 1:nrow(splitter)


## @knitr -----------break-------------


###
# Calculate correlations and kurtosis
###

## @knitr msit-dmn-correlation
sub_splitter <- subset(splitter, condition=="MSIT")
correlation_msit_dmn <- ddply(sub_splitter, .(subject, study, scan, run), function(sdf) {
    ts <- tss[[sdf$index]][,dmn]
    if (sdf$study == "CCB")
        r <- cor(ts, ccb_waver)
    else
        r <- cor(ts, ccd_waver)
    z <- atanh(r)
    c(msit.r=r, msit.z=z)    
})
# collapse across scan and run
correlation_msit_dmn <- ddply(correlation_msit_dmn, .(study, subject), numcolwise(mean))
# plot
ggplot(correlation_msit_dmn, aes(x=msit.z)) + 
    geom_histogram(aes(fill=..count..), binwidth=0.1) + 
    facet_grid(. ~ study) + 
    xlab("Correlations between MSIT task design and DN Signal (Fischer Z)")

## @knitr connectivity
connectivity_rest <- ddply(sub_splitter, .(study, subject, scan, run), function(sdf) {
    ts <- tss[[sdf$index]][,c(dmn,tps)]
    rs <- cor(ts)[-1,1]
    zs <- atanh(rs)
    data.frame(network=orig_network_names[tps], connectivity.r=rs, connectivity.z=zs)
})
# collapse across scan and run
connectivity_rest <- ddply(connectivity_rest, .(study, subject, network), numcolwise(mean))
# plot
ggplot(connectivity_rest, aes(x=connectivity.z)) + 
    geom_histogram(aes(fill=..count..), binwidth=0.1) + 
    facet_grid(network ~ study) + 
    xlab("Connectivity with DMN (Fischer Z)")

## @knitr combine
df <- merge(correlation_msit_dmn, connectivity_rest, by=c("study", "subject"))


## @knitr -----------break-------------


###
# Plot and significance
###

## @knitr functions
to_outlier <- function(x) factor((x>0.1)*1, levels=c(0,1), labels=c("yes", "no"))
wrap_lmrob <- function(f, df) {
    reg <- summary(lmrob(f, df, maxit.scale=500))
    print(reg)
    df$outlier <- to_outlier(reg$weights)
    df$weights <- reg$weights
    df$id <- 1:nrow(df)
    df
}

## @knitr signif
tmpdf <- ddply(df, .(network), function(sdf) {
    cat("\nConnectivity between DN and", as.character(sdf$network[1]), "\n")  
    wrap_lmrob(msit.z ~ connectivity.z, sdf)
})

## @knitr plots
grid <- ddply(df, .(network), function(sdf) {
    grid <- data.frame(
        connectivity.z = seq(min(sdf$connectivity.z), max(sdf$connectivity.z), length=20)
    )
    grid$msit.z <- predict(lmrob(msit.z ~ connectivity.z, tmpdf, maxit.scale=500), newdata=grid)
    grid
})
p0 <- ggplot(tmpdf, aes(x=connectivity.z, y=msit.z)) +
        geom_hline(aes(yintercept=0)) + 
        geom_vline(aes(yintercept=0), linetype='dashed') + 
        facet_grid(. ~ network) + 
        xlab("Connectivity with DN (Fischer Z)") + 
        ylab("Correlation between DN Signal\n& MSIT Task Design (Fischer Z)")
if (any(df$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(shape=1, size=8, color=brewer.pal(3,"Pastel1")[2]) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue")
} else {
    p <- p0 + 
            geom_point(shape=1, size=8, color=brewer.pal(3,"Pastel1")[2]) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
print(p)

