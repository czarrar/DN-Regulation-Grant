#!/usr/bin/env python

from pandas import *
import numpy as np
from os import path

behavdir = "/home2/data/Projects/CCD/behavior"

df = read_csv(path.join(behavdir, "ccd_totals.csv"), index_col=0)

# remove first 2 subjects since there is no imaging data on them
df = df[2:]
df.index = range(len(df))

# exclude subjects with BDI > 16 (i.e., those with moderate to high depression)
excluded_subjects = df[df['BDI']>16]['study_id']
excluded_subjects.tofile(path.join(behavdir, "excluded_subjects.txt"), sep="\n")

df = df[df['BDI']<=16]
df.index = range(0,len(df))
df.to_csv(path.join(behavdir, "ccd_totals_touse.csv"))
df.to_csv(path.join(behavdir, "ccd_totals_touse.xlsx"))
