library(plyr)
library(e1071)
library(ggplot2)

basedir <- "/home2/data/Projects/CCD"

network_names <- c("medial visual", "occipital pole visual", "lateral visual", "default network", "cerebellum", "sensorimotor", "auditory", "executive control", "right frontoparietal", "left frontoparietal")
nn <- gsub(" ", ".", network_names)
important_networks <- c(4,8,9,10)

###
# DMN-Kurtosis with MSIT-DMN
###

# Read in the phenotype data
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)[1:9,]

# Read in Task Waver
waver <- read.table(file.path(basedir, "rois/waver_msit_design.1D"))
waver <- as.vector(as.matrix(waver))

# Read in the RSNs
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
rest_tcs <- llply(fnames[1:9], read.table, .progress="text")
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_01/rsn10.1D")))
msit_run1_tcs <- llply(fnames, read.table, .progress="text")
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/msit/run_02/rsn10.1D")))
msit_run2_tcs <- llply(fnames, read.table, .progress="text")

# Correlate Waver with RSN
msit_run1_cors <- laply(msit_run1_tcs, cor, waver, .progress="text")
msit_run2_cors <- laply(msit_run2_tcs, cor, waver, .progress="text")
msit_cors <- (msit_run1_cors + msit_run2_cors)/2

# Kurtosis
rest_kurtosis <- laply(rest_tcs, apply, 2, kurtosis, .progress="text")
msit_run1_kurtosis <- laply(msit_run1_tcs, apply, 2, kurtosis, .progress="text")
msit_run2_kurtosis <- laply(msit_run2_tcs, apply, 2, kurtosis, .progress="text")
msit_kurtosis <- (msit_run1_kurtosis + msit_run2_kurtosis)/2

# Combine
df <- data.frame(
    id = rep(phenos$study_id, 2), 
    scan = factor(rep(c("rest", "msit"), each=9)), 
    kurtosis = c(rest_kurtosis[,4], msit_kurtosis[,4]), 
    msit_cors = rep(msit_cors[,4], 2)
)

# Plot
ggplot(df, aes(x=kurtosis, y=msit_cors)) + 
    geom_point() + 
    geom_smooth(method="lm") + 
    facet_grid(. ~ scan) + 
    labs(title="DMN-Kurtosis vs MSIT-DMN", x="Kurtosis", y="MSIT Task/BOLD Correlation")
ggsave(file.path(basedir, "analysis/combined/dmn-kurtosis_vs_msit-dmn.png"))

# Stats
summary(aov(msit_cors[,4] ~ rest_kurtosis[,4]))
## p = 0.4
summary(aov(msit_cors[,4] ~ msit_kurtosis[,4]))
## p = 0.07


###
# DMN-Kurtosis vs DMN-connectivity
###

# Read in the phenotype data
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names=1)

# Read in the RSNs
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
rest_tcs <- llply(fnames, read.table, .progress="text")

# Kurtosis
rest_kurtosis <- laply(rest_tcs, apply, 2, kurtosis, .progress="text")[,4]

# DMN Connectivity
rest_conn <- laply(rest_tcs, cor, .progress="text")[,,4]

# Combine
df <- data.frame(
    id = rep(phenos$study_id, 3), 
    connection = rep(c("Executive Control", "Right FrontoParietal", "Left FrontoParietal"), each=nrow(phenos)), 
    kurtosis = rep(rest_kurtosis, 3), 
    connectivity = c(rest_conn[,8:10])
)

# Plot
ggplot(df, aes(x=kurtosis, y=connectivity)) + 
    geom_point() + 
    geom_smooth(method="lm") + 
    facet_grid(. ~ connection) + 
    labs(title="DMN-Kurtosis vs DMN-Connectivity", 
            x="Kurtosis", y="Connectivity with DMN")
ggsave(file.path(basedir, "analysis/combined/dmn-kurtosis_vs_dmn-connectivity.png"))

# Stats
summary(aov(rest_conn[,8] ~ rest_kurtosis))
summary(aov(rest_conn[,9] ~ rest_kurtosis))
summary(aov(rest_conn[,10] ~ rest_kurtosis))
## all non-significant


###
# MSIT-DMN vs DMN-connectivity
###

# phenotype
phenos <- phenos[1:9,]

# msit-dmn
msit_cors <- msit_cors[,4]

# dmn-connectivity
rest_conn <- rest_conn[1:9,8:10]

# Combine
df <- data.frame(
    id = rep(phenos$study_id, 3), 
    connection = rep(c("Executive Control", "Right FrontoParietal", "Left FrontoParietal"), each=nrow(phenos)), 
    msit_cor = rep(msit_cors, 3), 
    connectivity = c(rest_conn)
)

# Plot
ggplot(df, aes(x=msit_cor, y=connectivity)) + 
    geom_point() + 
    geom_smooth(method="lm") + 
    facet_grid(. ~ connection) + 
    labs(title="MSIT-DMN vs DMN-Connectivity", 
            x="MSIT Task/Bold Correlation", y="Connectivity with DMN")
ggsave(file.path(basedir, "analysis/combined/msit-dmnn_vs_dmn-connectivity.png"))

# Stats
summary(aov(rest_conn[,1] ~ msit_cors))
summary(aov(rest_conn[,2] ~ msit_cors))
summary(aov(rest_conn[,2] ~ msit_cors))
## all non-significant
