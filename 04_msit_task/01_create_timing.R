#!/usr/env/bin Rscript --vanilla

# CCB with 224vols (TR=1.75s) structure is:
## 15 rest TRs
## 8 blocks (alternating coherent and incoherent) with 24 TRs each
## 17 rest TRs

begin_rest <- rep(0,15)
end_rest <- rep(0,17)

coherent <- c(begin_rest, rep(rep(c(1,0), each=24), 4), end_rest)
write.table(coherent, row.names=F, col.names=F, quote=F, file="CCB_coherent.1D")

incoherent <- c(begin_rest, rep(rep(c(0,1), each=24), 4), end_rest)
write.table(incoherent, row.names=F, col.names=F, quote=F, file="CCB_incoherent.1D")


# CCD with 153vols (TR=2s) structure is:
## 13 rest TRs
## 6 blocks (alternating coherent and incoherent) with 21 TRs each
## 14 rest TRs

begin_rest <- rep(0,13)
end_rest <- rep(0,14)

coherent <- c(begin_rest, rep(rep(c(1,0), each=21), 3), end_rest)

incoherent <- c(begin_rest, rep(rep(c(0,1), each=21), 3), end_rest)


