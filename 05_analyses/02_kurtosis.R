## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
library(robustbase)
basedir <- dirname(dirname(getwd())) # assume running in current directory
scriptdir <- dirname(getwd())

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")

## @knitr phenotypes
fname <- file.path(scriptdir, "data/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)

## @knitr timeseries
load(file.path(scriptdir, "data/ccb+ccd_time_series.rda"))
splitter <- attr(tss, 'split_labels')
splitter$index <- 1:nrow(splitter)
# Rest
sub_splitter <- subset(splitter, condition=="REST" & study=="CCD")
sub_splitter$subject <- factor(sub_splitter$subject)
rest_tcs <- daply(sub_splitter, .(subject), function(sdf) {
    tcs <- sapply(tss[sdf$index], function(x) x[,dmn])
    tc <- rowMeans(tcs)
    tc
})
# MSIT Run 1
sub_splitter <- subset(splitter, condition=="MSIT" & study=="CCD" & run==1)
sub_splitter$subject <- factor(sub_splitter$subject)
msit_run1_tcs <- daply(sub_splitter, .(subject), function(sdf) {
    tcs <- sapply(tss[sdf$index], function(x) x[,dmn])
    tc <- rowMeans(tcs)
    tc
})
# MSIT Run 2
sub_splitter <- subset(splitter, condition=="MSIT" & study=="CCD" & run==2)
sub_splitter$subject <- factor(sub_splitter$subject)
msit_run2_tcs <- daply(sub_splitter, .(subject), function(sdf) {
    tcs <- sapply(tss[sdf$index], function(x) x[,dmn])
    tc <- rowMeans(tcs)
    tc
})

## @knitr kurtosis
rest_kurtosis <- aaply(rest_tcs, 1, kurtosis, .progress="text")
msit_run1_kurtosis <- aaply(msit_run1_tcs, 1, kurtosis, .progress="text")
msit_run2_kurtosis <- aaply(msit_run2_tcs, 1, kurtosis, .progress="text")


## @knitr combine
# Rest
df.rest <- data.frame(
    phenos, 
    kurtosis=rest_kurtosis
)
# MSIT
n <- length(msit_run1_kurtosis)
df.msit <- data.frame(
    study_id = factor(rep(phenos$study_id[1:n], 2)), 
    run = factor(rep(c(1,2), each=n)), 
    rbind(phenos[1:n,-1], phenos[1:n,-1]), 
    kurtosis=c(msit_run1_kurtosis, msit_run2_kurtosis)
)
# MSIT (average 2 runs) & Rest
df.msit_and_rest <- data.frame(
    study_id = factor(rep(phenos$study_id[1:n], 2)), 
    scan = factor(rep(c("msit", "rest"), each=n)), 
    rbind(phenos[1:n,-1], phenos[1:n,-1]), 
    kurtosis=c((msit_run1_kurtosis+msit_run2_kurtosis)/2, rest_kurtosis[1:n])
)
mat.msit_and_rest <- data.frame(
    study_id = phenos$study_id[1:n], 
    msit = (msit_run1_kurtosis+msit_run2_kurtosis)/2, 
    rest = rest_kurtosis[1:n]
)
# Difference between MSIT (average 2 runs) & Rest
df.msit_vs_rest <- data.frame(
    study_id = phenos$study_id[1:n], 
    phenos[1:n,-1], 
    kurtosis=(msit_run1_kurtosis+msit_run2_kurtosis)/2 - rest_kurtosis[1:n]
)

## @knitr dist
# Rest
ggplot(df.rest, aes(x=kurtosis)) + 
    geom_histogram(binwidth=0.2) + 
    labs(title="Default-Mode Network During Rest", x="Kurtosis")
# MSIT
ggplot(df.msit, aes(x=kurtosis)) + 
    geom_histogram(binwidth=0.2) + 
    labs(title="Default-Mode Network During MSIT", x="Kurtosis")
# Scan Effect (MSIT/Rest)
ggplot(mat.msit_and_rest, aes(x=msit, y=rest)) + 
    geom_point() + 
    geom_smooth(method="lm") + 
    labs(title="Default-Mode Network Kurtosis", x="MSIT", y="Rest")
# MSIT vs Rest
ggplot(df.msit_vs_rest, aes(x=kurtosis)) + 
    geom_histogram(binwidth=0.2) + 
    labs(title="Default-Mode Network MSIT vs Rest", x="Kurtosis")


## @knitr anova-totals
# During Rest
summary(lmrob(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.rest, maxit.scale=500))
# During MSIT
summary(lmrob(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.msit, maxit.scale=500))
# Effects of Scan (MSIT/Rest)
summary(lmrob(kurtosis ~ Age + Sex + scan + Error(study_id), df.msit_and_rest, maxit.scale=500))
# MSIT vs REST
summary(lmrob(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.msit_vs_rest, maxit.scale=500))

## @knitr anova-subscales
# During Rest
summary(lmrob(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.rest, maxit.scale=500))
# During MSIT
summary(lmrob(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.msit, maxit.scale=500))
# MSIT vs REST
summary(lmrob(kurtosis ~ Age + Sex + SIPI + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.msit_vs_rest, maxit.scale=500))

## @knitr bimodal
library(mixtools)

mixmdl <- normalmixEM(rest_kurtosis, k=2)
plot(mixmdl, which=2)

grps <- factor((mixmdl$posterior[,1]>0.5)*1+1, labels=c("Positive", "Negative"))
df.pos <- df.rest[grps=="Positive",]
df.neg <- df.rest[grps=="Negative",]

msg = paste(
    "Mean Kurtosis of:\n", 
    "positive group = ", mean(df.pos$kurtosis), "\n", 
    "negative group = ", mean(df.neg$kurtosis), "\n"
)
cat(msg)

## @knitr bi-anova-totals
# Negative Group
summary(lmrob(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.neg, maxit.scale=500))
# Positive Group
summary(lmrob(kurtosis ~ Age + Sex + SIPI, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + RRS, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + ERQ, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + BDI, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + AIM, df.pos, maxit.scale=500))

## @knitr bi-anova-subscales
# Negative Group
summary(lmrob(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.neg, maxit.scale=500))
# Positive Group
summary(lmrob(kurtosis ~ Age + Sex + RRS_Brooding, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + RRS_Reflection, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + RRS_Depression, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + PANAS_Positive, df.pos, maxit.scale=500))
summary(lmrob(kurtosis ~ Age + Sex + PANAS_Negative, df.pos, maxit.scale=500))
