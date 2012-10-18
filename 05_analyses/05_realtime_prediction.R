## @knitr general-setup
library(plyr)
library(reshape)
library(e1071)
library(ggplot2)
library(vegan)
library(bcp)
library(RColorBrewer)
library(robustbase)
library(MASS)
basedir <- "/home2/data/Projects/CCD"
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10


## @knitr -----------break-------------

## @knitr phenotypes
fname <- file.path(basedir, "behavior/ccd_totals.csv")
phenos <- read.csv(fname, row.names=1)
phenos <- phenos[14:27,][-8,]  # CCD014 ... CCD027

## @knitr connectivity
# Read in time-series
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
tcs <- tcs[12:24,,] # CCD014 ... CCD027
# Calculate correlations
rest_conn_all <- aaply(tcs, 1, cor)
rest_conn <- rest_conn_all[,tps,1]   # only look at DMN connectivity with TP networks
colnames(rest_conn) <- network_names[tps]
names(dimnames(rest_conn)) <- c("subjects", "networks")
# Mean
colMeans(rest_conn)

## @knitr kurtosis
# only for DMN
rest_kurtosis <- aaply(tcs[,,4], 1, kurtosis)
# distribution
ggplot(data.frame(x=rest_kurtosis), aes(x=x)) +
    geom_histogram(binwidth=0.2) +
    geom_hline(aes(yintercept=0)) + 
    labs(x="Kurtosis")

## @knitr autocorrelation
# only for DMN
rest_lags <- aaply(tcs[,,4], 1, function(vec) {
    acors <- acf(vec, plot=F)
    lags <- c(acors$lag)
    acors <- c(acors$acf)
    for (i in 1:length(lags)) {
        if (acors[i] < 0)
            break
    }
    lags[i]
})
# distribution
ggplot(data.frame(x=rest_lags), aes(x=x)) +
    geom_histogram(binwidth=5) +
    geom_hline(aes(yintercept=0)) + 
    labs(x="Number of Lags for Zero Autocorrelation")

## @knitr changepoints
rest_changes <- aaply(tcs[,,4], 1, function(vec) {
    bcp.0 <- bcp(vec)
    sum(bcp.0$posterior.prob>0.5, na.rm=T)
})
# distribution
ggplot(data.frame(x=rest_changes), aes(x=x)) +
    geom_histogram(binwidth=5) +
    geom_hline(aes(yintercept=0)) + 
    labs(x="Kurtosis")
# sample subject
plot(bcp(tcs[1,,4]))

## @knitr prediction
fname <- file.path(basedir, "behavior/CCD_full_dframe.csv")
preds <- read.csv(fname)
preds <- preds[,c(2,6,10)] # only look at DMN
colnames(preds)[3] <- "R"
preds$Z <- atanh(preds$R)
# Plot distribution
ggplot(preds, aes(x=R)) +
    geom_histogram(binwidth=0.1) + 
    geom_vline(aes(xintercept=0), linetype='dashed') + 
    geom_hline(aes(yintercept=0)) + 
    facet_grid(ScanType ~ .) + 
    labs(title="Real-Time Prediction Accuracies", x="Correlation")


## @knitr -----------break-------------

## @knitr prediction-difference
# Plot
mat <- cast(preds, Subject ~ ScanType, value="R")
mat$diff <- apply(mat[,2:3], 1, diff)
ggplot(mat, aes(x=diff)) +
    geom_histogram(binwidth=0.1) + 
    geom_vline(aes(xintercept=0), linetype='dashed') + 
    geom_hline(aes(yintercept=0)) + 
    labs(title="Real-Time Prediction Accuracies", x="Difference between FB and NoFB")
# Significance
t.test(Z ~ ScanType, preds, paired=T)

## @knitr combine
# Moving forward only look at Feedback condition
preds <- subset(preds, ScanType=="FB")
# Combine phenotype info, kurtosis, and prediction accuracy
df <- data.frame(
    phenos, 
    prediction = preds$Z, 
    kurtosis = rest_kurtosis,  
    lag = rest_lags, 
    nchanges = rest_changes, 
    rest_conn
)
colnames(df)[[1]] <- "Subject"
# Needed functions
to_outlier <- function(x) factor((x>0.1)*1, levels=c(0,1), labels=c("yes", "no"))
wrap_lmrob <- function(f, df) {
    reg <- summary(lmrob(f, df, maxit.scale=500))
    print(reg)
    df$outlier <- to_outlier(reg$weights)
    df$weights <- reg$weights
    df
}

## @knitr kurtosis-prediction
# Significance
tdf <- wrap_lmrob(prediction ~ kurtosis, df)
# Plot
p <- ggplot(tdf, aes(x=kurtosis, y=prediction)) +
        geom_vline(aes(xintercept=0), linetype='dashed') + 
        geom_hline(aes(yintercept=0)) + 
        labs(x="DMN Kurtosis", y="DN Regulation (Fischer Z)")
if (any(tdf$outlier=="yes")) {
    p <- p + geom_point(aes(color=outlier)) + 
            geom_smooth(data=tdf[tdf$outlier=="no",], method="lm") + 
            scale_color_manual(values=c("black","red"))
} else {
    p <- p + geom_point() + geom_smooth(method="lm")
}
print(p)

## @knitr lag-prediction
tdf <- wrap_lmrob(prediction ~ lag, df)
# Plot
p <- ggplot(tdf, aes(x=lag, y=prediction)) +
        geom_vline(aes(xintercept=0), linetype='dashed') + 
        geom_hline(aes(yintercept=0)) + 
        xlab("Number of Lags for Zero Autocorrelation") + 
        ylab("DN Regulation (Fischer Z)")
if (any(tdf$outlier=="yes")) {
    p <- p + geom_point(aes(color=outlier)) + 
            geom_smooth(data=tdf[tdf$outlier=="no",], method="lm") + 
            scale_color_manual(values=c("black","red"))
} else {
    p <- p + geom_point() + geom_smooth(method="lm")
}
print(p)

## @knitr nchanges-prediction
tdf <- wrap_lmrob(prediction ~ nchanges, df)
# Plot
p <- ggplot(tdf, aes(x=nchanges, y=prediction)) +
        geom_vline(aes(xintercept=0), linetype='dashed') + 
        geom_hline(aes(yintercept=0)) + 
        xlab("Number of Change Points in DMN Time-Series") + 
        ylab("DN Regulation (Fischer Z)")
if (any(tdf$outlier=="yes")) {
    p <- p + geom_point(aes(color=outlier)) + 
            geom_smooth(data=tdf[tdf$outlier=="no",], method="lm") + 
            scale_color_manual(values=c("black","red"))
} else {
    p <- p + geom_point() + geom_smooth(method="lm")
}
print(p)

## @knitr connectivity-prediction
# Combine
tmpdf <- data.frame(
    df[rep(1:nrow(df), length(network_names[tps])),c("Subject", "Age", "Sex", "prediction")], 
    network = rep(network_names[tps], each=nrow(df)), 
    connectivity = atanh(c(rest_conn))
)
# Outliers
tmpdf <- ddply(tmpdf, .(network), function(sdf) {
    wrap_lmrob(connectivity ~ Age + Sex + prediction, sdf)
})
# Plot
ggplot(tmpdf, aes(x=connectivity, y=prediction, shape=network)) + 
    geom_vline(aes(xintercept=0), linetype='dashed') + 
    geom_hline(aes(yintercept=0)) + 
    geom_point() + 
    geom_smooth(method="lm") + 
    facet_grid(network ~ .) + 
    labs(x="Connectivity with DMN (Fischer Z)", y="DN Regulation (Fischer Z)")

## @knitr cwas-prediction
# CWAS between DMN and all networks
d <- as.dist(1-cor(t(rest_conn_all[,-4,4])))
adonis(d ~ Age + Sex + prediction, df, permutations=4999)
# CWAS between DMN and TP networks
d <- as.dist(1-cor(t(rest_conn)))
adonis(d ~ Age + Sex + prediction, df, permutations=4999)

## @knitr functions-prediction
brainbehavior <- function(names) {
    # Significance
    f <- paste("prediction ~ Age + Sex +", paste(names, collapse=" + "))
    f <- as.formula(f)
    tdf <- wrap_lmrob(f, df)
    
    # Reorganize
    tdf$id <- 1:nrow(tdf)
    bb.df <- ddply(tdf, .(Subject), function(sdf) {
        sdf <- data.frame(
            sdf[rep(1,length(names)), c("id", "Subject","prediction","outlier","weights")], 
            measure = names, 
            behavior = as.numeric(sdf[,names])
        )
        sdf
    })
    
    ## Get best fit line
    #model <- lmrob(prediction ~ measure, bb.df)
    #grid <- ddply(bb.df, .(measure), function(sdf) {
    #    data.frame(
    #        behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20), 
    #        measure=rep(sdf$measure[1], 20)
    #    )
    #})
    #grid$prediction <- predict(model, newdata=grid)
    
    # Plot
    p <- ggplot(bb.df, aes(x=behavior, y=prediction, label=id)) +
            geom_hline(aes(yintercept=0)) + 
            ylim(0,0.6) + 
            xlab("Scale Score") + 
            ylab("DN Regulation (Fischer Z)") + 
            facet_grid(. ~ measure, scales="free_x")
    if (any(bb.df$outlier=="yes")) {
        p <- p + 
                geom_point(data=bb.df[bb.df$outlier=="yes",], size=8, 
                            color=brewer.pal(3,"Pastel1")[1]) +
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(size=5) +
                geom_smooth(method=rlm) + 
                scale_color_discrete(name="Measure")
    } else {
        p <- p + 
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(size=5) +
                geom_smooth(method=rlm) + 
                scale_color_discrete(name="Measure")
    }
    p
}

## @knitr totals-prediction
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior(names)

## @knitr sipi-prediction
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior(names)

## @knitr erq-prediction
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior(names)

## @knitr rrs-prediction
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior(names)

## @knitr panas-prediction
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior(names)


