#!/usr/env/bin Rscript --vanilla

# CCD with 153vols structure is:
## 13 rest TRs
## 6 blocks (alternating coherent and incoherent) with 21 TRs each
## 14 rest TRs

begin_rest <- rep(0,13)
end_rest <- rep(0,14)

coherent <- c(begin_rest, rep(rep(c(1,0), each=21), 3), end_rest)

incoherent <- c(begin_rest, rep(rep(c(0,1), each=21), 3), end_rest)


# CCD with 224vols structure is:
## 15 rest TRs
## 6 blocks (alternating coherent and incoherent) with 32 TRs each
## 17 rest TRs

begin_rest <- rep(0,15)
end_rest <- rep(0,17)

coherent <- c(begin_rest, rep(rep(c(1,0), each=32), 3), end_rest)
write.table(coherent, row.names=F, col.names=F, quote=F, file="CCB_coherent.1D")

incoherent <- c(begin_rest, rep(rep(c(0,1), each=32), 3), end_rest)
write.table(incoherent, row.names=F, col.names=F, quote=F, file="CCB_incoherent.1D")
