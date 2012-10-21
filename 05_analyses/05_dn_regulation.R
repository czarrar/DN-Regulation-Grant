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
basedir <- dirname(dirname(getwd()))
datadir <- file.path(basedir, "scripts/data")
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10


## @knitr -----------break-------------

## @knitr phenotypes
fname <- file.path(datadir, "ccd_totals.csv")
phenos <- read.csv(fname, row.names=1)
phenos <- phenos[14:27,][-c(8,13),]  # CCD014 ... CCD027 (NO CCD021 and CCD026)

## @knitr connectivity
# Read in time-series
load(file.path(datadir, "ccb+ccd_time_series_all.rda"))
splitter <- attr(tss, 'split_labels')
splitter$index <- 1:nrow(splitter)
sub_splitter <- subset(splitter, subject %in% toupper(phenos$study_id))
tcs <- laply(tss[sub_splitter$index], function(x) x)
# Calculate correlations
rest_conn_all <- aaply(tcs, 1, cor)
rest_conn <- rest_conn_all[,tps,dmn]   # only look at DMN connectivity with TP networks
colnames(rest_conn) <- network_names[tps]
names(dimnames(rest_conn)) <- c("subjects", "networks")
# Mean
colMeans(rest_conn)

## @knitr kurtosis
# only for DMN
rest_kurtosis <- aaply(tcs[,,dmn], 1, kurtosis)
# distribution
ggplot(data.frame(x=rest_kurtosis), aes(x=x)) +
    geom_histogram(binwidth=0.2) +
    geom_hline(aes(yintercept=0)) + 
    labs(x="Kurtosis")

## @knitr autocorrelation
# only for DMN
rest_lags <- aaply(tcs[,,dmn], 1, function(vec) {
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
rest_changes <- aaply(tcs[,,dmn], 1, function(vec) {
    bcp.0 <- bcp(vec)
    sum(bcp.0$posterior.prob>0.5, na.rm=T)
})
# distribution
ggplot(data.frame(x=rest_changes), aes(x=x)) +
    geom_histogram(binwidth=5) +
    geom_hline(aes(yintercept=0)) + 
    labs(x="Kurtosis")
# sample subject
plot(bcp(tcs[1,,dmn]))

## @knitr prediction
fname <- file.path(datadir, "CCD_full_dframe.csv")
preds <- read.csv(fname)
preds <- preds[-c(23,24),c(2,6,10)] # no CCD026; only look at DMN
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
# Save
write.csv(df, file="z_data1.csv")
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
            geom_smooth(data=tdf[tdf$outlier=="no",], method=rlm) + 
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
            geom_smooth(data=tdf[tdf$outlier=="no",], method=rlm) + 
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
            geom_smooth(data=tdf[tdf$outlier=="no",], method=rlm) + 
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
    wrap_lmrob(prediction ~ Age + Sex + connectivity, sdf)
})
# Plot
ggplot(tmpdf, aes(x=connectivity, y=prediction, shape=network)) + 
    geom_vline(aes(xintercept=0), linetype='dashed') + 
    geom_hline(aes(yintercept=0)) + 
    geom_point() + 
    geom_smooth(method=rlm) + 
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
brainbehavior.multiple <- function(names, with_age_sex=TRUE) {
    # Significance
    if (with_age_sex)
        f <- paste("prediction ~ Age + Sex +", paste(names, collapse=" + "))
    else
        f <- paste("prediction ~", paste(names, collapse=" + "))
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
    
    # Get best fit line
    model <- lmrob(prediction ~ behavior + measure, bb.df, maxit.scale=500)
    grid <- ddply(bb.df[bb.df$outlier=="no",], .(measure), function(sdf) {
        data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20), 
            measure=rep(sdf$measure[1], 20)
        )
    })
    grid$prediction <- predict(model, newdata=grid)
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=prediction)) +
            geom_hline(aes(yintercept=0)) + 
            ylim(0,0.6) + 
            xlab("Scale Score") + 
            ylab("DN Regulation (Fischer Z)") + 
            facet_grid(. ~ measure, scales="free_x")
    if (any(bb.df$outlier=="yes")) {
        p <- p0 + 
                geom_point(data=bb.df[bb.df$outlier=="yes",], size=8, 
                            color=brewer.pal(3,"Pastel1")[1]) +
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(aes(label=id), size=5) +
                geom_line(data=grid, color="blue") + 
                scale_color_discrete(name="Measure")
    } else {
        p <- p0 + 
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(aes(label=id), size=5) +
                geom_line(data=grid, color="blue") + 
                scale_color_discrete(name="Measure")
    }
    p
}
brainbehavior.single <- function(names, with_age_sex=TRUE) {
    # Significance
    bb.df <- ldply(names, function(name) {
        cat("\nRunning regression for", name, "\n")
        if (with_age_sex)
            f <- paste("prediction ~ Age + Sex +", name)
        else
            f <- paste("prediction ~", name)
        f <- as.formula(f)
        tdf <- wrap_lmrob(f, df)
        tdf$id <- 1:nrow(tdf)
        tdf$measure <- name
        tdf$behavior <- tdf[[name]]
        cat("\n")
        tdf[,c("id", "Subject", "measure", "behavior", "prediction", "outlier", "weights")]
    })
    bb.df$measure <- factor(bb.df$measure)
    bb.df$outlier <- factor(bb.df$outlier)
    
        
    # Get best fit line
    grid <- ddply(bb.df[bb.df$outlier=="no",], .(measure), function(sdf) {
        model <- lmrob(prediction ~ behavior, sdf, maxit.scale=500)
        sgrid <- data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20)
        )
        sgrid$prediction <- predict(model, newdata=sgrid)
        sgrid$measure <- sdf$measure[1]
        sgrid
    })
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=prediction)) +
            geom_hline(aes(yintercept=0)) + 
            ylim(0,0.6) + 
            xlab("Scale Score") + 
            ylab("DN Regulation (Fischer Z)") + 
            facet_grid(. ~ measure, scales="free_x")
    if (any(bb.df$outlier=="yes")) {
        p <- p0 + 
                geom_point(data=bb.df[bb.df$outlier=="yes",], size=8, 
                            color=brewer.pal(3,"Pastel1")[1]) +
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(aes(label=id), size=5) +
                geom_line(data=grid, color="blue") + 
                scale_color_discrete(name="Measure")
    } else {
        p <- p0 + 
                geom_point(aes(color=measure), shape=1, size=8) +
                geom_text(aes(label=id), size=5) +
                geom_line(data=grid, color="blue") + 
                scale_color_discrete(name="Measure")
    }
    p
}

## @knitr multiple-totals-prediction
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.multiple(names)

## @knitr multiple-totals-prediction-no-bdi
names <- c("SIPI", "RRS", "ERQ", "AIM", "PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names)

## @knitr multiple-sipi-prediction
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.multiple(names)

## @knitr multiple-erq-prediction
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.multiple(names)

## @knitr multiple-rrs-prediction
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior.multiple(names)

## @knitr multiple-panas-prediction
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names)

## @knitr single-totals-prediction
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.single(names)

## @knitr single-sipi-prediction
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.single(names)

## @knitr single-erq-prediction
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.single(names)

## @knitr single-rrs-prediction
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
brainbehavior.single(names)

## @knitr single-panas-prediction
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.single(names)


## @knitr multiple-totals-prediction-just-say-no
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-totals-prediction-no-bdi-just-say-no
names <- c("SIPI", "RRS", "ERQ", "AIM", "PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-sipi-prediction-just-say-no
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-erq-prediction-just-say-no
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-rrs-prediction-just-say-no
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-panas-prediction-just-say-no
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, with_age_sex=FALSE)


## @knitr single-totals-prediction-just-say-no
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-sipi-prediction-just-say-no
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-erq-prediction-just-say-no
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-rrs-prediction-just-say-no
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-panas-prediction-just-say-no
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.single(names, with_age_sex=FALSE)


## @knitr multiple-totals-prediction-just-say-no-no1
df <- df[-1,]
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-totals-prediction-no-bdi-just-say-no-no1
names <- c("SIPI", "RRS", "ERQ", "AIM", "PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-sipi-prediction-just-say-no-no1
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-erq-prediction-just-say-no-no1
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-rrs-prediction-just-say-no-no1
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior.multiple(names, with_age_sex=FALSE)

## @knitr multiple-panas-prediction-just-say-no-no1
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, with_age_sex=FALSE)


## @knitr single-totals-prediction-just-say-no-no1
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-sipi-prediction-just-say-no-no1
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-erq-prediction-just-say-no-no1
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-rrs-prediction-just-say-no-no1
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
brainbehavior.single(names, with_age_sex=FALSE)

## @knitr single-panas-prediction-just-say-no-no1
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.single(names, with_age_sex=FALSE)
