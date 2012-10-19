#!/usr/bin/env Rscript

# This should be run after 03_setup_level2.py
# it deals with a rank-deficient model

df <- read.csv("z_predesign.csv")

df$age <- scale(df$age, scale=F)
df$study <- factor(df$study)
df$scan <- factor(df$scan)

formula <- . ~ age + sex + study + scan + subject
formula[[2]] <- NULL
X.frame <- model.frame(formula, df, drop.unused.levels = TRUE)
X <- model.matrix(formula, X.frame)
X <- X[,colSums(X)!=1]
qrhs <- qr(X)
if (qrhs$rank < ncol(X))
    cat("Warning: model is rank deficient, adjusting\n")
X <- X[, qrhs$pivot, drop = FALSE]
X <- X[, 1:qrhs$rank, drop = FALSE]

write.table(X, row.names=F, col.names=F, quote=F, file="z_design_forfsl.txt")
write.csv(X, row.names=F, quote=F, file="z_design.csv")
