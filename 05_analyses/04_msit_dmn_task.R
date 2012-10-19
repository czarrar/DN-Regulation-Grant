## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
library(RColorBrewer)
library(robustbase)
basedir <- "/home2/data/Projects/CCD"
scriptdir <- file.path(basedir, "scripts/04_msit_task")
oldtheme <- theme_set(theme_bw())

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- 4
tps <- 8:10


## @knitr -----------break-------------


###
# Setup
###

## @knitr subject-info
phenos <- read.csv(file.path(basedir, "scripts/data/ccd_totals_touse.csv"))[1:9,-1] # only CCD subjects with MSIT
subinfo <- read.csv(file.path(scriptdir, "z_predesign.csv"))
subinfo$study <- factor(subinfo$study, labels=c("CCB", "CCD"))
subinfo$scan <- factor(subinfo$scan)
subinfo$run <- factor(subinfo$run)
subinfo$sex <- factor(subinfo$sex)
# Remove CCD participants with high errors
subinfo <- subinfo[!((subinfo$subject=="CCD004"|subinfo$subject=="CCD008")&(subinfo$run==1)),]
# Plot Age
ggplot(subinfo, aes(x=age, fill=..count..)) + 
    geom_histogram(binwidth=5) + 
    facet_grid(study ~ .)
# Plot Sex
ggplot(subinfo, aes(x=sex, fill=study)) + 
    geom_bar()

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
# Group Average Plot
###

## @knitr msit-design
# CCB
ccb_coherent   <- read.table(file.path(scriptdir, "CCB_coherent.1D"))[,1]
ccb_incoherent <- read.table(file.path(scriptdir, "CCB_incoherent.1D"))[,1]
tmp <- diff(ccb_coherent - ccb_incoherent)
# +1 to adjust for relative; +3 to adjust for HRF delay
ccb_event_tpts <- data.frame(
    xmin = c(0, which(tmp!=0) + 1 + 3), 
    xmax = c(which(tmp!=0) + 1 + 3, length(ccb_coherent)), 
    block = factor(c("Fixation", rep(c("Coherent", "Incoherent"), length.out=sum(tmp!=0)-1), "Fixation"))
)
ccb_event_tpts$xmin <- ccb_event_tpts$xmin * 1.75
ccb_event_tpts$xmax <- ccb_event_tpts$xmax * 1.75
# CCD
ccd_coherent   <- read.table(file.path(scriptdir, "CCD_coherent.1D"))[,1]
ccd_incoherent <- read.table(file.path(scriptdir, "CCD_incoherent.1D"))[,1]
tmp <- diff(ccd_coherent - ccd_incoherent)
# +1 to adjust for relative; +3 to adjust for HRF delay
ccd_event_tpts <- data.frame(
    xmin = c(0, which(tmp!=0) + 1 + 3), 
    xmax = c(which(tmp!=0) + 1 + 3, length(ccd_coherent)), 
    block = factor(c("Fixation", rep(c("Coherent", "Incoherent"), length.out=sum(tmp!=0)-1), "Fixation"))
)
ccd_event_tpts$xmin <- ccd_event_tpts$xmin * 2
ccd_event_tpts$xmax <- ccd_event_tpts$xmax * 2

## @knitr msit-dmn-average
# CCB
sub_splitter <- subset(splitter, condition=="MSIT" & study=="CCB")
sub_splitter$subject <- factor(sub_splitter$subject)
ccb_msit_tcs <- daply(sub_splitter, .(subject), function(sdf) {
    tcs <- sapply(tss[sdf$index], function(x) x[,dmn])
    tc <- rowMeans(tcs)
    tc
})
ccb_msit_tc_ave <- data.frame(
    timepoint = (1:ncol(ccb_msit_tcs)) * 1.75, 
    bold = colMeans(ccb_msit_tcs)
)
# CCD
sub_splitter <- subset(splitter, condition=="MSIT" & study=="CCD")
sub_splitter$subject <- factor(sub_splitter$subject)
ccd_msit_tcs <- daply(sub_splitter, .(subject), function(sdf) {
    tcs <- sapply(tss[sdf$index], function(x) x[,dmn])
    tc <- rowMeans(tcs)
    tc
})
ccd_msit_tc_ave <- data.frame(
    timepoint = (1:ncol(ccd_msit_tcs)) * 2, 
    bold = colMeans(ccd_msit_tcs)
)

## @knitr msit-dmn-average-plot
# CCB
p <- ggplot(ccb_msit_tc_ave) +
        geom_rect(data=ccb_event_tpts, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=block)) +
        scale_fill_manual(name="", breaks=c("Fixation", "Coherent", "Incoherent"), 
                            values=brewer.pal(4, "Pastel2")[-4]) + 
        geom_line(aes(x=timepoint, y=bold), color="darkblue", size=0.75) + 
        scale_x_discrete(name="Time (secs)", limits=c(0,224*1.75), expand=c(0,0)) + 
        scale_y_continuous(name="BOLD Signal", limits=c(-0.3,0.3), expand=c(0,0))
print(p)      
# CCD
p <- ggplot(ccd_msit_tc_ave) +
        geom_rect(data=ccd_event_tpts, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=block)) +
        scale_fill_manual(name="", breaks=c("Fixation", "Coherent", "Incoherent"), 
                            values=brewer.pal(4, "Pastel2")[-4]) + 
        geom_line(aes(x=timepoint, y=bold), color="darkblue", size=0.75) + 
        scale_x_discrete(name="Time (secs)", limits=c(0,153*2), expand=c(0,0)) + 
        scale_y_continuous(name="BOLD Signal", limits=c(-0.3,0.3), expand=c(0,0))
print(p)      


## @knitr -----------break-------------


###
# Task-DMN Correlation
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


## @knitr -----------break-------------


###
# Task-DMN Correlation with Phenotypes (only with CCD)
###

## @knitr combine
names(phenos)[1:3] <- c("subject", "age", "sex")
phenos$subject <- factor(toupper(phenos$subject))
phenos$sex <- factor(phenos$sex, labels=c("F", "M"))
df <- merge(correlation_msit_dmn, phenos)
names(df)[[3]] <- "correlation.dmn_with_task"

## @knitr functions-phenotype
to_outlier <- function(x) factor((x>0.1)*1, levels=c(0,1), labels=c("yes", "no"))
wrap_lmrob <- function(f, df) {
    reg <- summary(lmrob(f, df, maxit.scale=500))
    print(reg)
    df$outlier <- to_outlier(reg$weights)
    df$weights <- reg$weights
    df
}
brainbehavior.multiple <- function(names, df) {
    # Significance
    f <- paste("correlation.dmn_with_task ~ age + sex +", paste(names, collapse=" + "))
    f <- as.formula(f)
    tdf <- wrap_lmrob(f, df)
    
    # Reorganize
    tdf$id <- 1:nrow(tdf)
    bb.df <- ddply(tdf, .(subject), function(sdf) {
        sdf <- data.frame(
            sdf[rep(1,length(names)), c("id", "subject","correlation.dmn_with_task","outlier","weights")], 
            measure = names, 
            behavior = as.numeric(sdf[,names])
        )
        sdf
    })
    
    # Get best fit line
    model <- lmrob(correlation.dmn_with_task ~ behavior + measure, bb.df, maxit.scale=500)
    grid <- ddply(bb.df, .(measure), function(sdf) {
        data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20), 
            measure=rep(sdf$measure[1], 20)
        )
    })
    grid$correlation.dmn_with_task <- predict(model, newdata=grid)
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=correlation.dmn_with_task)) +
            geom_hline(aes(yintercept=0)) + 
            xlab("Scale Score") + 
            ylab("Correlation between DMN Signal\n& Task Design (Fischer Z)") + 
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
brainbehavior.single <- function(names, df) {
    # Significance
    bb.df <- ldply(names, function(name) {
        cat("\nRunning regression for", name, "\n")
        f <- paste("correlation.dmn_with_task ~ age + sex +", name)
        f <- as.formula(f)
        tdf <- wrap_lmrob(f, df)
        tdf$id <- 1:nrow(tdf)
        tdf$measure <- name
        tdf$behavior <- tdf[[name]]
        cat("\n")
        tdf[,c("id", "subject", "measure", "behavior", "correlation.dmn_with_task", "outlier", "weights")]
    })
    bb.df$measure <- factor(bb.df$measure)
    bb.df$outlier <- factor(bb.df$outlier)
    
        
    # Get best fit line
    grid <- ddply(bb.df, .(measure), function(sdf) {
        model <- lmrob(correlation.dmn_with_task ~ behavior, sdf, maxit.scale=500)
        sgrid <- data.frame(
            behavior=seq(min(sdf$behavior), max(sdf$behavior), length=20)
        )
        sgrid$correlation.dmn_with_task <- predict(model, newdata=sgrid)
        sgrid$measure <- sdf$measure[1]
        sgrid
    })
    
    # Plot
    p0 <- ggplot(bb.df, aes(x=behavior, y=correlation.dmn_with_task)) +
            geom_hline(aes(yintercept=0)) + 
            xlab("Scale Score") + 
            ylab("Correlation between DMN Signal\n& Task Design (Fischer Z)") + 
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

## @knitr multiple-totals-correlation-bdi
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.multiple(names, df)

## @knitr multiple-totals-correlation-panas
names <- c("SIPI", "RRS", "ERQ", "AIM", "PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, df)

## @knitr multiple-sipi-correlation
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.multiple(names, df)

## @knitr multiple-erq-correlation
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.multiple(names, df)

## @knitr multiple-rrs-correlation
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior.multiple(names, df)

## @knitr multiple-panas-correlation
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.multiple(names, df)

## @knitr single-totals-correlation
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.single(names, df)

## @knitr single-sipi-correlation
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.single(names, df)

## @knitr single-erq-correlation
names <- c("ERQ_Reappraisal", "ERQ_Suppression")
brainbehavior.single(names, df)

## @knitr single-rrs-correlation
names <- c("RRS_Brooding", "RRS_Depression", "RRS_Reflection")
brainbehavior.single(names, df)

## @knitr single-panas-correlation
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.single(names, df)


