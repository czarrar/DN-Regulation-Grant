## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
basedir <- "/home2/data/Projects/CCD"

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")

## @knitr phenotypes
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)

## @knitr timeseries

# Rest
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
rest_tcs <- laply(fnames, function(f) as.matrix(read.table(f)), .progress="text")

# MSIT Run 1
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_01/rsn10.1D")))
msit_run1_tcs <- laply(fnames, function(f) as.matrix(read.table(f)), .progress="text")

# MSIT Run 2
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_02/rsn10.1D")))
msit_run2_tcs <- laply(fnames, function(f) as.matrix(read.table(f)), .progress="text")

## @knitr only-dmn
rest_tcs <- rest_tcs[,,dmn]
msit_run1_tcs <- msit_run1_tcs[,,dmn]
msit_run2_tcs <- msit_run2_tcs[,,dmn]

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
summary(aov(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.rest))
# During MSIT
summary(aov(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.msit))
# Effects of Scan (MSIT/Rest)
summary(aov(kurtosis ~ Age + Sex + scan + Error(study_id), df.msit_and_rest))
# MSIT vs REST
summary(aov(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.msit_vs_rest))

## @knitr anova-subscales
# During Rest
summary(aov(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.rest))
# During MSIT
summary(aov(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.msit))
# MSIT vs REST
summary(aov(kurtosis ~ Age + Sex + SIPI + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.msit_vs_rest))

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
summary(aov(kurtosis ~ Age + Sex + SIPI + RRS + ERQ + BDI + AIM, df.neg))
# Positive Group
summary(aov(kurtosis ~ Age + Sex + SIPI, df.pos))
summary(aov(kurtosis ~ Age + Sex + RRS, df.pos))
summary(aov(kurtosis ~ Age + Sex + ERQ, df.pos))
summary(aov(kurtosis ~ Age + Sex + BDI, df.pos))
summary(aov(kurtosis ~ Age + Sex + AIM, df.pos))

## @knitr bi-anova-subscales
# Negative Group
summary(aov(kurtosis ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative, df.neg))
# Positive Group
summary(aov(kurtosis ~ Age + Sex + RRS_Brooding, df.pos))
summary(aov(kurtosis ~ Age + Sex + RRS_Reflection, df.pos))
summary(aov(kurtosis ~ Age + Sex + RRS_Depression, df.pos))
summary(aov(kurtosis ~ Age + Sex + PANAS_Positive, df.pos))
summary(aov(kurtosis ~ Age + Sex + PANAS_Negative, df.pos))
