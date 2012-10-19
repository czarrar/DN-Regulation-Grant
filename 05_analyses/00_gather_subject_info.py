#!/usr/bin/env python

from glob import glob
from os import path
from pandas import *
import numpy as np

basepath = "/home/data/Projects/CCD"


###
# CCB
###

analysis = path.join(basepath, "analysis/subjects")
ccb_filt = read_csv(path.join(basepath, "scripts/data/ccb_subject_info_filt.csv"))
ccb_filt['study'] = "CCB"
ccb_filt = ccb_filt[["subject", "age", "sex", "study", "scan", "condition", "run"]]


###
# CCD
###

subjects = np.array([ ['CCD%03i' % i]*3 for i in xrange(3, 12) ]).flatten()
runs = np.array([ [1, 2, 1] for i in xrange(3,12) ]).flatten()
conditions = np.array([ ["MSIT", "MSIT", "REST"] for i in xrange(3,12) ]).flatten()
ccd_all = DataFrame({'subject': subjects, 'condition': conditions, 'run': runs})
ccd_all['scan'] = 1
ccd_all['study'] = "CCD"

# exclude CCD004 run 1; CCD008 run 1
bad_inds = ((ccd_all.subject == "CCB004") | (ccd_all.subject == "CCB008")) & (ccd_all.run == 1) & (ccd_all.condition == "MSIT")
ccd_filt0 = ccd_all[~bad_inds]

# get demographics
ccd_filt1 = read_csv(path.join(basepath, "behavior/ccd_all_filtered.csv"))[["study_id", "Age", "Sex"]]
ccd_filt1.columns = ["subject", "age", "sex"]
ccd_filt1.subject = [ s.upper() for s in ccd_filt1.subject ]

# merge demographics with run info and save
ccd_filt = merge(ccd_filt0, ccd_filt1, on="subject", how="inner")
ccd_filt = ccd_filt[["subject", "age", "sex", "study", "scan", "condition", "run"]]
ccd_filt.sex[ccd_filt.sex=="Male"] = "M"
ccd_filt.sex[ccd_filt.sex=="Female"] = "F"


###
# Combined
###

filt = ccb_filt.copy()
filt = filt.append(ccd_filt, ignore_index=True)
filt.to_csv(path.join(basepath, "scripts/data/ccb+ccd_filtered.csv"), index=False)
