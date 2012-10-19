#!/usr/bin/env python

# Select the CCB subjects to include in any analyses

from os import path
from pandas import *
import numpy as np
import re

basepath = "/home2/data/Projects/CCD"
inpath = path.join(basepath, "CCB/extras")
outpath = path.join(basepath, "scripts/data")

###
# Read in demographics
###
print "Reading in demographics"
demo = read_csv(path.join(inpath, "CCB_demo.csv"))
demo_df = demo.ix[:,0:4]
demo_df.columns = ["subject", "age", "sex", "scanner"]


###
# Read in motion info for Rest and MSIT
###
print "Reading in motion info"
motion = read_csv(path.join(inpath, "CCB_motion_corrected.csv"))
# Parser for filenames
parse_motion_fname = re.compile("rp_(?P<subject>CCB[0-9]{3})_scan(?P<scan>[12])_(?P<condition>[A-Z]+)(?P<run>[12]).*.1D")
# Organize info
cols = ['subject', 'condition', 'scan', 'run', 'max_motion']
d = { x : [] for x in cols }
for i,row in motion.T.iteritems():
    details = parse_motion_fname.search(row['File']).groupdict()
    d['subject'].append(details['subject'])
    d['condition'].append(details['condition'])
    d['scan'].append(int(details['scan']))
    d['run'].append(int(details['run']))
    d['max_motion'].append(float(row['Max Motion (mm)']))
motion_df = DataFrame(d, columns=cols)


###
# Read in errors for MSIT
###
print "Reading in errors for MSIT"
errors = read_table(path.join(inpath, "CCB_MSIT_Errors.txt"), header=None, sep=' ')
errors.columns = ["filename", "errors", "noresponses"]
# Parser for filenames
parse_behav_fname = re.compile("(?P<subject>CCB[0-9]{3})_(?:Scan)?(?P<scan>2)?[_]?msit_(?P<run>[12])_behavior.1D")
# Organize info
cols = ['subject', 'condition', 'scan', 'run', 'errors']
d = { x : [] for x in cols}
for i,row in errors.T.iteritems():
    details = parse_behav_fname.search(row['filename']).groupdict()
    if details['scan'] is None: details['scan'] = '1'
    d['subject'].append(details['subject'])
    d['condition'].append('MSIT')
    d['scan'].append(int(details['scan']))
    d['run'].append(int(details['run']))
    d['errors'].append(int(row['errors']))
errors_df = DataFrame(d, columns=cols)


###
# Join together the different tables
###
print "Joining the 3 tables together and saving"
left = merge(demo_df, motion_df, on='subject', how='outer')
df_all = merge(left, errors_df, on=['subject', 'condition', 'scan', 'run'], how='outer')
df_all.to_csv(path.join(outpath, "ccb_subject_info_all.csv"), index=False)


###
# Determine what stays and what goes
# - motion must be less than 3
# - errors must be less than 10
###
print "Excluding subjects with motion <3 and errors >10 and saving"
inds = (df_all.max_motion<3) & ((df_all.errors<10) | (df_all.condition=="REST"))
df_filt = df_all[inds]
df_filt.to_csv(path.join(outpath, "ccb_subject_info_filt.csv"), index=False)
df_bad = df_all[~inds]
df_bad.to_csv(path.join(outpath, "ccb_subject_info_bad.csv"), index=False)
