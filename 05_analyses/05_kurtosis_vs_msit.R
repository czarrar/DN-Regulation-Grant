## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
library(robustbase)
library(RColorBrewer)
basedir <- "/home2/data/Projects/CCD"
scriptdir <- file.path(basedir, "scripts/04_msit_task")
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10


## @knitr -----------break-------------


###
# Setup
###

## @knitr subject-info
fname <- file.path(scriptdir, "z_predesign.csv")
subinfo <- read.csv(fname)
subinfo$study <- factor(subinfo$study, labels=c("CCB", "CCD"))
subinfo$sex <- factor(subinfo$sex)
# Remove CCD participants with high errors
subinfo <- subinfo[!((subinfo$subject=="CCD004"|subinfo$subject=="CCD008")&(subinfo$run==1)),]
# Plot Age
ggplot(subinfo, aes(x=age, fill=..count..)) + geom_histogram(binwidth=5) + facet_grid(study ~ .)
# Plot Sex
ggplot(subinfo, aes(x=sex, fill=study)) + geom_bar()

## @knitr task-waver
# CCB
ccb_waver <- read.table(file.path(scriptdir, "level1_ccb_template.mat"), skip=5)
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
ccd_waver <- read.table(file.path(scriptdir, "level1_ccd_template.mat"), skip=5)
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
load(file.path(basedir, "scripts/data/ccb+ccd_time_series.rda"))
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
    c(r=r, z=z)    
})
# collapse across scan and run
correlation_msit_dmn <- ddply(correlation_msit_dmn, .(subject), numcolwise(mean))

## @knitr kurtosis
sub_splitter <- subset(splitter, condition=="REST")
kurtosis_rest <- ddply(sub_splitter, .(subject, study, scan, run), function(sdf) {
    ts <- tss[[sdf$index]][,dmn]
    c(kurtosis=kurtosis(ts))
})
# collapse across scan and run
kurtosis_rest <- ddply(kurtosis_rest, .(subject), numcolwise(mean))

## @knitr combine
df <- merge(correlation_msit_dmn, kurtosis_rest, by='subject')
# plot msit-dmn
ggplot(df, aes(x=z, fill=..count..)) + 
    geom_histogram(binwidth=5) + 
    labs(x="Correlation between MSIT task design & DN signal", y="Count")
# plot kurtosis
ggplot(df, aes(x=kurtosis, fill=..count..)) + 
    geom_histogram(binwidth=5) + 
    labs(x="Kurtosis During Rest in DN", y="Count")



## @knitr -----------break-------------


###
# Plot and significance
###

## @knitr comparison-functions
to_outlier <- function(x) factor((x>0.1)*1, levels=c(0,1), labels=c("yes", "no"))
wrap_lmrob <- function(f, df) {
    reg <- summary(lmrob(f, df, maxit.scale=500))
    print(reg)
    df$outlier <- to_outlier(reg$weights)
    df$weights <- reg$weights
    df$id <- 1:nrow(df)
    df
}

## @knitr comparison-signif
tmpdf <- wrap_lmrob(z ~ kurtosis, df)

## @knitr comparison-plot
grid <- data.frame(
    kurtosis = seq(min(df$kurtosis), max(df$kurtosis), length=20)
)
grid$z <- predict(lmrob(z ~ kurtosis, tmpdf, maxit.scale=500), newdata=grid)
p0 <- ggplot(tmpdf, aes(x=kurtosis, y=z)) +
        geom_hline(aes(yintercept=0)) + 
        geom_vline(aes(yintercept=0), linetype='dashed') + 
        xlab("DN Kurtosis During Rest") + 
        ylab("Correlation between DN Signal\n& MSIT Task Design")
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
p
