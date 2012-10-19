#!/usr/bin/env python

from glob import glob
from os import path
from pandas import *
import numpy as np

basepath = "/home/data/Projects/CCD"


###
# CCB: Generate paths to MSIT feats
###

analysis = path.join(basepath, "analysis/subjects")
ccb_filt = read_csv(path.join(basepath, "scripts/data/ccb_subject_info_filt.csv"))
ccb_filt = ccb_filt[ccb_filt.condition=="MSIT"]
ccb_filt['study'] = 1
ccb_filt = ccb_filt[["subject", "age", "sex", "study", "scan", "run"]]

def gen_path(row):
    find_path = path.join(analysis, row['subject'], 
        "msit_%02i" % row['scan'], 
        "run_%02i" % row['run'], 
        "task_hrf.feat"
    )
    paths = glob(find_path)
    if len(paths) == 0:
        print "Error: no paths found for %s" % find_path
        print row
        return None
    if len(paths) > 1:
        print "Error: %i paths found for %s" % (len(paths), find_path)
        print row
    return paths[0]

ccb_filepaths = ccb_filt.apply(gen_path, axis=1)


###
# CCD: Generate paths to MSIT feats
###

analysis = path.join(basepath, "analysis/subjects")

# setup data frame
subjects = np.array([ ['CCD%03i' % i]*2 for i in xrange(3, 12) ]).flatten()
runs = np.array([ [1, 2] for i in xrange(3,12) ]).flatten()
ccd_all = DataFrame({'subject': subjects, 'run': runs}, columns=['subject', 'run'])
ccd_all['scan'] = 1

# exclude CCD004 run 1; CCD008 run 1
bad_inds = ((ccd_all.subject == "CCB004") & (ccd_all.run == 1)) | ((ccd_all.subject == "CCB008") & (ccd_all.run == 1))
ccd_filt0 = ccd_all[~bad_inds]

# get demographics
ccd_filt1 = read_csv(path.join(basepath, "behavior/ccd_all_filtered.csv"))[["study_id", "Age", "Sex"]]
ccd_filt1.columns = ["subject", "age", "sex"]
ccd_filt1.subject = [ s.upper() for s in ccd_filt1.subject ]

# merge demographics with run info and save
ccd_filt = merge(ccd_filt0, ccd_filt1, on="subject", how="inner")
ccd_filt['study'] = 2
ccd_filt = ccd_filt[["subject", "age", "sex", "study", "scan", "run"]]
ccd_filt.sex[ccd_filt.sex=="Male"] = "M"
ccd_filt.sex[ccd_filt.sex=="Female"] = "F"
ccd_filt.to_csv(path.join(basepath, "scripts/data/ccd_subject_info_filt.csv"), index=False)

# now generate paths
def gen_path(row):
    find_path = path.join(analysis, 
        row['subject'], "msit", 
        "run_%02i" % row['run'], 
        "task_hrf.feat"
    )
    paths = glob(find_path)
    if len(paths) == 0:
        print "Error: no paths found for %s" % find_path
        print row
        return None
    if len(paths) > 1:
        print "Error: %i paths found for %s" % (len(paths), find_path)
        print row
    return paths[0]

ccd_filepaths = ccd_filt.apply(gen_path, axis=1)


###
# Combine CCB + CCD
###

filt = ccb_filt.copy()
filt = filt.append(ccd_filt, ignore_index=True)
filt.to_csv(path.join("z_predesign.csv"), index=False)

filepaths = Series(ccb_filepaths.tolist() + ccd_filepaths.tolist())
filepaths.to_csv("z_funcpaths.txt", sep=" ", index=False)


###
# Create design file
###

design              = DataFrame(index=range(filt.shape[0]))
design['intercept'] = [1]*filt.shape[0]
design['age']       = (filt.age - filt.age.mean()).tolist()
design['sex']       = ((filt.sex == "M")*1).tolist()
design['study']     = (filt.study - 1).tolist()
design['scan']      = (filt.scan - 1).tolist()
runs = DataFrame(np.array([ ((filt.subject==name)*1).tolist() for name,group in filt.groupby('subject') ]).T, 
                    columns=[ "sub%02i" % (i+1) for i in range(len(filt.subject.unique())) ])
runs = runs.ix[:,0:-1] # exclude last subject since have intercept
design[runs.columns] = runs
design.to_csv("z_design.csv", index=False)
design.to_csv("z_design_forfsl.txt", header=False, index=False, sep=" ")

