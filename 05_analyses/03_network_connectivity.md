# Network Connectivity

Note that actual code is loaded from a different file.


```r
read_chunk("03_network_connectivity.R")
```



## Setup


```r
library(plyr)
library(reshape)
```

```
## Attaching package: 'reshape'
```

```
## The following object(s) are masked from 'package:plyr':
## 
## rename, round_any
```

```r
library(ggplot2)
```

```
## Loading required package: methods
```

```r
basedir <- "/home2/data/Projects/CCD"
```



```r
network_names <- c("medial visual", "occipital pole visual", "lateral visual", 
    "default network", "cerebellum", "sensorimotor", "auditory", "executive control", 
    "right frontoparietal", "left frontoparietal")
network_names <- gsub(" ", ".", network_names)
dmn <- which(network_names == "default.network")
tps <- 8:10
```


### Read in Data


```r
fname <- file.path(basedir, "behavior/ccd_totals_touse.csv")
phenos <- read.csv(fname, row.names = 1)
```



```r
fnames <- sort(Sys.glob(file.path(basedir, "analysis/subjects/*/rest/run_01/rsn10.1D")))
rest_tcs <- laply(fnames, function(f) as.matrix(read.table(f)))
rest_tcs <- rest_tcs[, , c(dmn, tps)]
```


Note since analyses will focus on the default-mode network and the task positive networks, time-series will be restricted to only those networks.


## Connectivity

For each subject, calculate the correlation between the DMN time-series and each of 3 TP time-series.


```r
rest_conn <- aaply(rest_tcs, 1, cor)
rest_conn <- rest_conn[, -1, 1]  # only look at DMN connectivity with TP networks
colnames(rest_conn) <- network_names[tps]
names(dimnames(rest_conn)) <- c("subjects", "networks")
rest_conn
```

```
##         networks
## subjects executive.control right.frontoparietal left.frontoparietal
##       1           0.083891            -0.191988             0.01978
##       2           0.084717             0.158241             0.29943
##       3           0.481685            -0.333323            -0.05599
##       4           0.473743             0.224075            -0.23917
##       5           0.455081             0.309532             0.25369
##       6           0.181667            -0.131148            -0.31457
##       7           0.268529            -0.034411             0.01492
##       8           0.294297            -0.039687            -0.18061
##       9           0.391924             0.201970             0.47545
##       10          0.424079            -0.010988             0.11693
##       11          0.121752            -0.030745             0.03444
##       12          0.040938             0.048803            -0.28633
##       13          0.422123             0.280951             0.02877
##       14          0.182321             0.141492            -0.13037
##       15          0.219550             0.072720            -0.03550
##       16          0.146454            -0.008729             0.03150
##       17         -0.008234            -0.066889             0.15641
##       18          0.078349             0.173039             0.24086
##       19          0.218871            -0.154438            -0.15270
##       20          0.220138             0.199491            -0.29388
##       21          0.390617             0.028109            -0.08993
##       22          0.175987            -0.113726            -0.07558
##       23          0.162978            -0.284650            -0.23095
```


### Distribution


```r
df <- melt(rest_conn)
colnames(df) <- c("subject", "network", "connectivity")
ggplot(df, aes(x = connectivity)) + geom_histogram(binwidth = 0.1) + facet_grid(network ~ 
    .) + geom_vline(aes(xintercept = 0), linetype = "dashed") + geom_hline(aes(yintercept = 0)) + 
    labs(title = "Connectivity with DMN", x = "Correlation")
```

![plot of chunk connectivity-distribution](figure/connectivity-distribution.png) 



## Brain-Behavior

Generally the approach below is to examine the relationship between each phenotypic measure and DMN/TP network connectivity via scatter plots and a table of Z-Scores assessing the significance of the relationships.

### Needed Functions











