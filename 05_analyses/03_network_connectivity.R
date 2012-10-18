## @knitr general-setup
library(plyr)
library(reshape)
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
phenos <- read.csv(fname, row.names=1)

## @knitr timeseries
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
rest_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
rest_tcs <- rest_tcs[,,c(dmn,tps)]


## @knitr -----------break-------------


## @knitr connectivity
rest_conn <- aaply(rest_tcs, 1, cor)
rest_conn <- rest_conn[,-1,1] # only look at DMN connectivity with TP networks
colnames(rest_conn) <- network_names[tps]
names(dimnames(rest_conn)) <- c("subjects", "networks")
rest_conn

## @knitr connectivity-distribution
df <- melt(rest_conn)
colnames(df) <- c("subject", "network", "connectivity")
ggplot(df, aes(x=connectivity)) + 
    geom_histogram(binwidth=0.1) + 
    facet_grid(network ~ .) + 
    geom_vline(aes(xintercept=0), linetype="dashed") +
    geom_hline(aes(yintercept=0)) + 
    labs(title="Connectivity with DMN", x="Correlation")


## @knitr -----------break-------------

## @knitr brain-behavior-functions
# note: i am being lazy here and using global variables inside these function calls
brainbehavior.plot <- function(names) {
    zscore <- function(x) (x-mean(x))/sd(x)
    bb.df <- ddply(df, .(subject, network), function(sdf) {
        sid <- sdf$subject[[1]]
        sdf <- data.frame(
            sdf[rep(1,length(names)),], 
            measure = names, 
            behavior = as.numeric(phenos[sid,names])
        )
        sdf
    })
    bb.df <- ddply(bb.df, .(network, measure), function(sdf) {
        sdf$behavior <- zscore(sdf$behavior)
        sdf
    })
    
    p <- ggplot(bb.df, aes(x=behavior, y=connectivity)) + 
        geom_point(aes(color=measure)) +
        geom_smooth(aes(color=measure), method="lm") + 
        facet_grid(network ~ measure) + 
        labs(title="Connectivity-Phenotype Relationships", 
                x="Scale Z-Score", y="Correlation")
    print(p)
    
    invisible(bb.df)
}

brainbehavior.signif <- function(names) {
    brain.behavior <- aaply(rest_conn, 2, function(connectivity) {
        f <- paste("connectivity ~ Age + Sex +", paste(names, collapse=" + "))
        f <- as.formula(f)
        tab <- summary(aov(f, phenos))
        tab[[1]][names,"Pr(>F)"]
    })
    colnames(brain.behavior) <- names
    brain.behavior <- round(qt(brain.behavior, Inf, lower.tail=F), 2)
    brain.behavior
}

create_table <- function(tab) {
    tab <- cbind(rownames(tab), tab)
    cat(paste(colnames(tab), collapse=" | "), "\n")
    cat(paste(rep("---", ncol(tab)), collapse=" | "), "\n")
    cat(apply(tab, 1, function(x) paste(x, collapse=" | ")), sep = "\n")
}

## @knitr brain-behavior-totals
names <- c("SIPI", "RRS", "ERQ", "BDI", "AIM")
brainbehavior.plot(names)
create_table(brainbehavior.signif(names))
## TODO: indicate what is significant from above! (maybe in bold?)

## @knitr brain-behavior-sipi
names <- c("SIPI_PAC", "SIPI_GFFD", "SIPI_PCD")
brainbehavior.plot(names)
brainbehavior.signif(names)

## @knitr brain-behavior-rrs
names <- c("RRS_Brooding", "RRS_Reflection", "RRS_Depression")
brainbehavior.plot(names)
brainbehavior.signif(names)

## @knitr brain-behavior-panas
names <- c("PANAS_Positive", "PANAS_Negative")
brainbehavior.plot(names)
brainbehavior.signif(names)

