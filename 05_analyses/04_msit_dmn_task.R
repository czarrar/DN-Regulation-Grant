## @knitr general-setup
library(plyr)
library(e1071)
library(ggplot2)
basedir <- "/home2/data/Projects/CCD"

## @knitr network-setup
network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10


## @knitr -----------break-------------


## @knitr phenotypes
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)[1:9,] # only for 9 subjects with MSIT

# Read in Task Waver
task_waver <- read.table(file.path(basedir, "rois/waver_msit_design.1D"))
task_waver <- as.vector(as.matrix(task_waver))

## @knitr timeseries
# Run 1
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_01/rsn10.1D")))
msit_run1_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
msit_run1_tcs <- msit_run1_tcs[,,dmn]
# Run 2
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_02/rsn10.1D")))
msit_run2_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
msit_run2_tcs <- msit_run2_tcs[,,dmn]


## @knitr -----------break-------------


## @knitr correlate
msit_run1_cors <- aaply(msit_run1_tcs, 1, cor, task_waver)
msit_run2_cors <- aaply(msit_run2_tcs, 1, cor, task_waver)



BELOW NEEDS TO BE FINISHED


# Combine Stuff
df <- data.frame(
    id=factor(rep(phenos$study_id, 2)), 
    run=factor(rep(c("Run 1", "Run 2"), each=9)), 
    rbind(phenos[,-1], phenos[,-1]), 
    rbind(msit_run1_cors, msit_run2_cors)
)

# Run regression
res.totals <- laply(nn, function(network) {
    summary(aov(df[[network]] ~ Age + Sex + SIPI + RRS + ERQ + AIM + Error(id), df))
})
res.subtotals <- laply(nn, function(network) {
    summary(aov(df[[network]] ~ Age + Sex + RRS_Brooding + RRS_Reflection + RRS_Depression + PANAS_Positive + PANAS_Negative + Error(id), df))
})

# Save some of it

# Get average of signal in DMN
msit_run1_avetc <- rowMeans(sapply(msit_run1_tcs, function(x) x[,4]))
msit_run2_avetc <- rowMeans(sapply(msit_run2_tcs, function(x) x[,4]))
msit_avetc <- (msit_run1_avetc + msit_run2_avetc)/2

# Plot
ggplot(df, aes(x=default.network)) + geom_histogram(binwidth=0.1) + facet_grid(. ~ run) +  geom_vline(aes(xintercept=0), linetype="dashed", color="red") +  labs(title="Default-Mode Network", x="Correlation between Task Regressor and BOLD Signal")
ggsave(file.path(basedir, "analysis/msit/dmn_hist_cor.png"))
zscore <- function(x) (x-mean(x))/sd(x)
pdf <- data.frame(
    time=seq(0,length.out=length(waver)),
    value=c(zscore(msit_avetc), zscore(waver)), 
    type=factor(rep(c("BOLD Signal", "Task"), each=length(waver))) 
)
ggplot(pdf, aes(x=time, y=value, color=type)) + geom_hline(aes(yintercept=0), linetype="dashed") + geom_line() + labs(title="Default-Mode Network", x="Time (secs)", y="Z-Score")
ggsave(file.path(basedir, "analysis/msit/dmn_task_signal.png"))
