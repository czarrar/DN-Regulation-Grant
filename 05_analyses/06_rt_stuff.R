## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
library(RColorBrewer)
library(robustbase)
basedir <- dirname(dirname(getwd()))    # assume running in current direcotry
scriptdir <- file.path(basedir, "scripts")
datadir <- file.path(basedir, "scripts/data")
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
orig_network_names <- c("Medial Visual", "Occipital Pole Visual", "Lateral Visual", "Default Network", "Cerebellum", "Sensorimotor", "Auditory", "Executive Control", "Right Frontoparietal", "Left Frontoparietal")
network_names <- gsub(" ", ".", tolower(orig_network_names))
dmn <- 4
tps <- 8:10


## @knitr -----------break-------------


###
# Setup
###

## @knitr subject-info
# Phenotypic Info (only 9 subjects in CCD)
phenos <- read.csv(file.path(datadir, "ccd_totals_touse.csv"))[1:9,-1] 
colnames(phenos)[1:3] <- c("subject", "age", "sex")
phenos$subject <- toupper(phenos$subject)
# Basic Subject Info for all
subinfo <- read.csv(file.path(datadir, "ccb+ccd_filtered.csv"))
subinfo$scan <- factor(subinfo$scan)
subinfo$run <- factor(subinfo$run)

## @knitr rts
load(file.path(datadir, "ccb+ccd_rts.rda")) # rts
rts <- rts[rts$rt>200,] # remove trials with RT < 200ms
rts$scan <- factor(rts$scan)
rts$run <- factor(rts$run)

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

## @knitr timeseries
# this loads the 'tss' object with attr(tss, 'split_labels') to get how stuff should be organized
load(file.path(basedir, "scripts/data/ccb+ccd_time_series.rda"))
splitter <- attr(tss, 'split_labels')
splitter$index <- 1:nrow(splitter)


## @knitr -----------break-------------

###
# Calculate
###

## @knitr rts-summarize
rts_df <- ddply(rts, .(study, subject, condition), function(sdf) {
    c(mean.rt=mean(sdf$rt), cv.rt=sd(sdf$rt)/mean(sdf$rt))
})
rts_df
# Plot Mean RTs
ggplot(rts_df, aes(x=mean.rt)) + 
    geom_histogram(aes(fill=..count..)) + 
    facet_grid(condition ~ study) + 
    xlab("Mean Reaction Time (msecs)")
# Plot RT Coefficient of Variation
ggplot(rts_df, aes(x=cv.rt)) + 
    geom_histogram(aes(fill=..count..)) + 
    facet_grid(condition ~ study) + 
    xlab("Reaction Time Coefficient of Variation")

## @knitr kurtosis
sub_splitter <- subset(splitter, condition=="REST")
kurtosis_rest <- ddply(sub_splitter, .(study, subject, scan, run), function(sdf) {
    ts <- tss[[sdf$index]][,dmn]
    c(kurtosis=kurtosis(ts))
})
# collapse across scan and run
kurtosis_rest <- ddply(kurtosis_rest, .(study, subject), numcolwise(mean))
# plot
ggplot(kurtosis_rest, aes(x=kurtosis)) + 
    geom_histogram(aes(fill=..count..)) + 
    facet_grid(. ~ study) + 
    xlab("Kurtosis")

## @knitr connectivity
connectivity_rest <- ddply(sub_splitter, .(study, subject, scan, run), function(sdf) {
    ts <- tss[[sdf$index]][,c(dmn,tps)]
    rs <- cor(ts)[-1,1]
    zs <- atanh(rs)
    data.frame(network=orig_network_names[tps], r=rs, z=zs)
})
# collapse across scan and run
connectivity_rest <- ddply(connectivity_rest, .(study, subject, network), numcolwise(mean))
# plot
ggplot(connectivity_rest, aes(x=z)) + 
    geom_histogram(aes(fill=..count..)) + 
    facet_grid(network ~ study) + 
    xlab("Connectivity with DMN (Fischer Z)")

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
# plot CCB
ggplot(correlation_msit_dmn[1:28,], aes(x=z, fill=..count..)) + 
    geom_histogram(binwidth=0.05) + 
    geom_hline(aes(yintercept=0)) + 
    geom_vline(aes(xintercept=0), linetype="dashed") + 
    xlab("DMN Correlation with Task Design, Incoherent > Coherent (Fischer Z)") + 
    ylab("Number of Subjects") + 
    ggtitle("CCB Dataset")
# plot CCD
ggplot(correlation_msit_dmn[29:37,], aes(x=z, fill=..count..)) + 
    geom_histogram(binwidth=0.1) + 
    geom_hline(aes(yintercept=0)) + 
    geom_vline(aes(xintercept=0), linetype="dashed") + 
    xlab("DMN Correlation with Task Design, Incoherent > Coherent (Fischer Z)") +
    ylab("Number of Subjects") + 
    ggtitle("CCD Dataset")

## @knitr combine
df_phenos <- merge(rts_df, phenos, by='subject')
df_kurtosis <- merge(rts_df, kurtosis_rest, by=c('study', 'subject'))
df_connectivity <- merge(rts_df, connectivity_rest, by=c('study', 'subject'))
df_msit_dn <- merge(rts_df, correlation_msit_dmn, by=c('subject'))
df_kurtosis_phenos <- merge(subset(df_kurtosis, study=="CCD", select=-1), phenos, by='subject')
df_connectivity_phenos <- merge(subset(df_connectivity, study=="CCD", select=-1), phenos, by='subject')


## @knitr -----------break-------------

###
# Plot and Significance with Brain Measures
###

## @knitr fone
to_outlier <- function(x) factor((x>0.1)*1, levels=c(0,1), labels=c("yes", "no"))
wrap_lmrob <- function(f, df) {
    reg <- summary(lmrob(f, df, maxit.scale=500))
    print(reg)
    df$outlier <- to_outlier(reg$weights)
    df$weights <- reg$weights
    df
}
get_grid <- function(df, f, y, x) {
    model <- lmrob(f, df, maxit.scale=500)
    grid <- data.frame(row.names=1:20)
    grid[[x]] <- seq(min(df[[x]]), max(df[[x]]), length=20)
    grid[[y]] <- predict(model, newdata=grid)
    grid
}

## @knitr kurtosis-meanrt
tmpdf <- ddply(df_kurtosis, .(condition), function(sdf) {
    cat("\nCondition:", as.character(sdf$condition[1]), "\n")
    wrap_lmrob(kurtosis ~ mean.rt, sdf)
})
tmpdf$id <- rep(1:length(unique(df_kurtosis$subject)), 2)
grid <- ddply(df_kurtosis, .(condition), get_grid, kurtosis ~ mean.rt, "kurtosis", "mean.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=mean.rt, y=kurtosis)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("Mean MSIT RT (msecs)") + 
        ylab("DN Kurtosis at Rest") + 
        facet_grid(. ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p

## @knitr kurtosis-cvrt
tmpdf <- ddply(df_kurtosis, .(condition), function(sdf) {
    cat("\nCondition:", as.character(sdf$condition[1]), "\n")
    wrap_lmrob(kurtosis ~ cv.rt, sdf)
})
tmpdf$id <- rep(1:length(unique(df_kurtosis$subject)), 2)
grid <- ddply(df_kurtosis, .(condition), get_grid, kurtosis ~ cv.rt, "kurtosis", "cv.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=cv.rt, y=kurtosis)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("MSIT RT Coefficient of Variation") + 
        ylab("DN Kurtosis at Rest") + 
        facet_grid(. ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p

## @knitr connectivity-meanrt
tmpdf <- ddply(df_connectivity, .(network, condition), function(sdf) {
    cat("\nConnectivity with", as.character(sdf$network[1]), "during", as.character(sdf$condition[1]), "\n")
    tmpdf <- wrap_lmrob(z ~ mean.rt, sdf)
    tmpdf$id <- 1:nrow(tmpdf)
    tmpdf
})
grid <- ddply(df_connectivity, .(network, condition), get_grid, z ~ mean.rt, "z", "mean.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=mean.rt, y=z)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("Mean MSIT RT (msecs)") + 
        ylab("Connectivity with DN at Rest (Fischer Z)") + 
        facet_grid(network ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p


## @knitr connectivity-cvrt
tmpdf <- ddply(df_connectivity, .(network, condition), function(sdf) {
    cat("\nConnectivity with", as.character(sdf$network[1]), "during", as.character(sdf$condition[1]), "\n")
    tmpdf <- wrap_lmrob(z ~ cv.rt, sdf)
    tmpdf$id <- 1:nrow(tmpdf)
    tmpdf
})
grid <- ddply(df_connectivity, .(network, condition), get_grid, z ~ cv.rt, "z", "cv.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=cv.rt, y=z)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("MSIT RT Coefficient of Variation") + 
        ylab("Connectivity with DN at Rest (Fischer Z)") + 
        facet_grid(network ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p

## @knitr msit_dn-meanrt
tmpdf <- ddply(df_msit_dn, .(condition), function(sdf) {
    cat("\nCondition:", as.character(sdf$condition[1]), "\n")
    wrap_lmrob(z ~ mean.rt, sdf)
})
tmpdf$id <- rep(1:length(unique(df_msit_dn$subject)), 2)
grid <- ddply(df_msit_dn, .(condition), get_grid, z ~ mean.rt, "z", "mean.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=mean.rt, y=z)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("Mean MSIT RT (msecs)") + 
        ylab("Correlation between MSIT Task Design\n& DN Signal (Fischer Z)") + 
        facet_grid(. ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p

## @knitr msit_dn-cvrt
tmpdf <- ddply(df_msit_dn, .(condition), function(sdf) {
    cat("\nCondition:", as.character(sdf$condition[1]), "\n")
    wrap_lmrob(z ~ cv.rt, sdf)
})
tmpdf$id <- rep(1:length(unique(df_msit_dn$subject)), 2)
grid <- ddply(df_msit_dn, .(condition), get_grid, z ~ cv.rt, "z", "cv.rt")
# Plot
p0 <- ggplot(tmpdf, aes(x=cv.rt, y=z)) +
        geom_hline(aes(yintercept=0)) +  
        xlab("MSIT RT Coefficient of Variation") + 
        ylab("Correlation between MSIT Task Design\n& DN Signal (Fischer Z)") + 
        facet_grid(. ~ condition, scales="free_x")
if (any(tmpdf$outlier=="yes")) {
    p <- p0 + 
            geom_point(data=tmpdf[tmpdf$outlier=="yes",], size=8, 
                        color=brewer.pal(3,"Pastel1")[1]) +
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
} else {
    p <- p0 + 
            geom_point(aes(color=condition), shape=1, size=8) +
            geom_text(aes(label=id), size=5) +
            geom_line(data=grid, color="blue") + 
            scale_color_discrete(name="Measure")
}
p


## @knitr -----------break-------------


###
# Plot and Significance with Phenotypic Measures
###

## @knitr ftwo
meanrt.single <- function(df, names, title) {
    # Significance
    bb.df <- ldply(names, function(name) {
        cat("\nRunning regression for", name, "\n")
        f <- paste("mean.rt ~", name)
        f <- as.formula(f)
        tdf <- wrap_lmrob(f, df)
        tdf$id <- 1:nrow(tdf)
        tdf$measure <- name
        tdf$behavior <- tdf[[name]]
        cat("\n")
        tdf[,c("id", "subject", "measure", "behavior", "mean.rt", "outlier", "weights")]
    })
    bb.df$measure <- factor(bb.df$measure)
    bb.df$outlier <- factor(bb.df$outlier)
    
    # Get best fit line
    grid <- ddply(bb.df, .(measure), function(sdf) {
        model <- lmrob(mean.rt ~ behavior, sdf, maxit.scale=500)
        sgrid <- data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20)
        )
        sgrid$mean.rt <- predict(model, newdata=sgrid)
        sgrid$measure <- sdf$measure[1]
        sgrid
    })
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=mean.rt)) +
            xlab("Scale Score") + 
            ylab("Mean MSIT RT (secs)") + 
            facet_grid(. ~ measure, scales="free_x") + 
            ggtitle(title)
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
cvrt.single <- function(df, names, title) {
    # Significance
    bb.df <- ldply(names, function(name) {
        cat("\nRunning regression for", name, "\n")
        f <- paste("cv.rt ~", name)
        f <- as.formula(f)
        tdf <- wrap_lmrob(f, df)
        tdf$id <- 1:nrow(tdf)
        tdf$measure <- name
        tdf$behavior <- tdf[[name]]
        cat("\n")
        tdf[,c("id", "subject", "measure", "behavior", "cv.rt", "outlier", "weights")]
    })
    bb.df$measure <- factor(bb.df$measure)
    bb.df$outlier <- factor(bb.df$outlier)
    
    # Get best fit line
    grid <- ddply(bb.df, .(measure), function(sdf) {
        model <- lmrob(cv.rt ~ behavior, sdf, maxit.scale=500)
        sgrid <- data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20)
        )
        sgrid$cv.rt <- predict(model, newdata=sgrid)
        sgrid$measure <- sdf$measure[1]
        sgrid
    })
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=cv.rt)) +
            xlab("Scale Score") + 
            ylab("MSIT RT Coefficient of Variation") + 
            facet_grid(. ~ measure, scales="free_x") + 
            ggtitle(title)
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

## @knitr select-coherent
title <- "Coherent Trials"
df_phenos_coherent <- subset(df_phenos, condition=="Coherent")

## @knitr coherent-meanrt-totals
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
meanrt.single(df_phenos_coherent, names, title)

## @knitr coherent-meanrt-sipi
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
meanrt.single(df_phenos_coherent, names, title)

## @knitr coherent-meanrt-erq
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
meanrt.single(df_phenos_coherent, names, title)

## @knitr coherent-meanrt-rrs
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
meanrt.single(df_phenos_coherent, names, title)

## @knitr coherent-meanrt-panas
names <- c("PANAS_Positive", "PANAS_Negative")
meanrt.single(df_phenos_coherent, names, title)

## @knitr coherent-cvrt-totals
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
cvrt.single(df_phenos_coherent, names, title)

## @knitr coherent-cvrt-sipi
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
cvrt.single(df_phenos_coherent, names, title)

## @knitr coherent-cvrt-erq
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
cvrt.single(df_phenos_coherent, names, title)

## @knitr coherent-cvrt-rrs
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
cvrt.single(df_phenos_coherent, names, title)

## @knitr coherent-cvrt-panas
names <- c("PANAS_Positive", "PANAS_Negative")
cvrt.single(df_phenos_coherent, names, title)

## @knitr select-incoherent
title <- "Incoherent Trials"
df_phenos_incoherent <- subset(df_phenos, condition=="Incoherent")

## @knitr incoherent-meanrt-totals
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
meanrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-meanrt-sipi
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
meanrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-meanrt-erq
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
meanrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-meanrt-rrs
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
meanrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-meanrt-panas
names <- c("PANAS_Positive", "PANAS_Negative")
meanrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-cvrt-totals
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
cvrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-cvrt-sipi
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
cvrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-cvrt-erq
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
cvrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-cvrt-rrs
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
cvrt.single(df_phenos_incoherent, names, title)

## @knitr incoherent-cvrt-panas
names <- c("PANAS_Positive", "PANAS_Negative")
cvrt.single(df_phenos_incoherent, names, title)
