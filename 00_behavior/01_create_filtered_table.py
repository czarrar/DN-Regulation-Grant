#!/usr/bin/env python

from pandas import *
import numpy as np
from os import path
from progressbar import ProgressBar, Percentage, Bar, ETA, Counter
import re

behavdir = "/home2/data/Projects/CCD/behavior"


###
# Functions
###

def filter_list(pattern, list_of_text):
    """returns list with elements that have matching pattern"""
    return Series([ text for text in list_of_text if re.search(pattern, text) ])

def new_cols_starting_with(prefix, n):
    """returns a series of column names starting with prefix"""
    return Series([ "%s_Q%02i" % (prefix, i+1) for i in range(n) ])

def compute_subscale_scores(scale, sdf, subscales):
    tdf = DataFrame(index=range(len(sdf)))
    for subscale,questions in subscales.iteritems():
        tdf["%s_%s_Total" % (scale, subscale)] = sdf.ix[:,questions-1].apply(sum, 1)
    return tdf

def append_subscales(df, sdf):
    """for some reason concat/append doesn't work properly"""
    for i,col in sdf.iteritems():
        df[col.name] = col
    return df
    
def compute_total_score(sdf):
    return sdf.apply(sum, 1)


###
# Load Data and Sort
###

rawdf = read_csv(path.join(behavdir, "ccd_all_raw.csv"))
rawdf = rawdf.sort("study_id")
rawdf.index = range(len(rawdf))

###
# Basics: id, age, sex
###

df = rawdf[["study_id", "Age", "Sex"]].sort("study_id")


###
# EHI: Edinburg Handedness Inventory
###

print "\nEHI"

ehi_recoding = {"Always Left": -50, "Usually Left": -25, 
                "No Preference": 0, 
                "Usually Right": 25, "Always Right": 50}
 
print "...copying and renaming the columns to be EHI_Q01...EHI_Q07"
raw_cols = filter_list("EHI", rawdf.columns)
ehi_cols = new_cols_starting_with("EHI", len(raw_cols))
df[ehi_cols] = rawdf[raw_cols]

print "...changing the text responses to a hard solid number"
df[ehi_cols] = df[ehi_cols].applymap(lambda x: ehi_recoding[x])

print "...getting total"
df["EHI_Total"] = compute_total_score(df[ehi_cols])


###
# SIPI: Short Imaginal Processes Inventory (with 3 sub-scales)
###

print "\nSIPI"

# SIPI questions to be negatively-keyed
sipi_negatively_keyed = Series(sorted([4,17,19,26,38,30,44,1,8,13,16,23,35,42]))

# SIPI questions associated with each scale
sipi_scales = {
    # Positive-Constructive Daydreaming
    "PCD": Series(sorted([2,7,10,15,22,28,31,36,40,43,4,17,19,26,38])), 
    # Guilt and Fear of Failure Daydreaming
    "GFFD": Series(sorted([3,6,9,12,14,18,21,24,27,32,34,37,41,30,44])), 
    # Poor Attentional Control
    "PAC": Series(sorted([5,11,20,25,29,33,39,45,1,8,13,16,23,35,42]))
}

print "...copying and renaming the columns to be SIPI_Q01...SIPI_Q45"
raw_cols = filter_list("SIPI", rawdf.columns)
sipi_cols = new_cols_starting_with("SIPI", len(raw_cols))
df[sipi_cols] = rawdf[raw_cols]

print "...reversing negatively-keyed items"
subcols = sipi_cols[sipi_negatively_keyed-1].tolist() # not sure why i need tolist()
df[subcols] = 6 - df[subcols]

print "...getting sub-scales"
sdf = compute_subscale_scores("SIPI", df[sipi_cols], sipi_scales)
df = append_subscales(df, sdf)

print "...getting total"
df["SIPI_Total"] = compute_total_score(df[sipi_cols])


###
# RRS: Ruminative Responses Scale
###

print "\nRRS"

# questions to extract RRS from the larger RSQ
# note: other scales of RSQ are not reliable or good predictors of depressio
rrs_qs_in_rsq = Series(sorted([5, 6, 7, 8, 14, 15, 18, 19, 21, 22, 25, 28, 30, 27, 40, 42, 43, 44, 46, 53, 56, 61]))

# Each scale with the associated RRS questions
rrs_scales = {
    "Brooding": Series([ 5, 10, 13, 15, 16]), # bad: 13, 
    "Depression": Series([ 1,  2,  3,  4,  6,  8,  9, 14, 17, 18, 19, 22]), 
    "Reflection": Series([ 7, 11, 12, 20, 21])
}

rrs_recoding = {"Almost Never": 1, "Sometimes": 2, 
                "Often": 3, "Almost Always": 4}

print "...copying and renaming the columns to be RRS_Q01...RRS_Q22"
raw_cols = filter_list("RSQ", rawdf.columns)[rrs_qs_in_rsq-1]
rrs_cols = new_cols_starting_with("RRS", len(raw_cols))
df[rrs_cols] = rawdf[raw_cols.tolist()] # not sure why i need the tolist

print "...changing the text responses to numbers"
df[rrs_cols] = df[rrs_cols].applymap(lambda x: rrs_recoding[x])

print "...getting sub-scales"
sdf = compute_subscale_scores("RRS", df[rrs_cols], rrs_scales)
df = append_subscales(df, sdf)

print "...getting total"
df["RRS_Total"] = compute_total_score(df[rrs_cols])


###
# ERQ: Emotion Regulation Questionnaire
###

print "\nERQ"

# Each scale with the associated ERQ questions
erq_scales = {
    "Reappraisal": Series([1, 3, 5, 7, 8, 10]), 
    "Suppression": Series([2, 4, 6, 9])
}

print "...copying and renaming the columns to be ERQ_Q01...ERQ_Q10"
raw_cols = filter_list("ERQ", rawdf.columns)
erq_cols = new_cols_starting_with("ERQ", len(raw_cols))
df[erq_cols] = rawdf[raw_cols]

print "...cleaning up the item scores (ugh)"
extract_number = lambda score: int(re.search("[1-9]", score).group())
df[erq_cols] = df[erq_cols].applymap(extract_number)

print "...getting sub-scales"
sdf = compute_subscale_scores("ERQ", df[erq_cols], erq_scales)
df = append_subscales(df, sdf)

print "...getting total"
df["ERQ_Total"] = df[erq_cols].apply(sum, 1)


###
# BDI: Beck Depression Inventory
###

print "\nBDI"

nquestions = 21

print "...copying and renaming the columns to be BDI_Q01...BDI_Q21"
raw_cols = filter_list("BDI", rawdf.columns)
bdi_cols = new_cols_starting_with("BDI", nquestions)

print "...complicated scoring"
def get_item_score(answers, questions):
    # Which answer did the person select
    choice = np.argwhere(answers).flatten()[-1] # take highest numbered choice
    # What is the score for each answer
    scores = np.zeros(len(questions))
    for i,question in enumerate(questions):
        search = re.search("\[([1-9]).*\]$", question)
        if search is None:
            scores[i] = 0
        else:
            scores[i] = search.group(1)
    # Based on person's choice what is the score for this item
    return scores[choice]
for i,name in enumerate(bdi_cols):
    question_cols = filter_list("%02i" % (i+1), raw_cols)
    df[name] = rawdf[question_cols].apply(get_item_score, axis=1, questions=question_cols)

print "...getting total"
df["BDI_Total"] = df[bdi_cols].apply(sum, 1)


###
# PANAS: Positive Affective and Negative Affective Schedule
###

print "\nPANAS"

# Each scale with the associated PANAS questions
panas_scales = {
    "Positive": Series([1, 3, 5, 9, 10, 12, 14, 16, 17, 19]), 
    "Negative": Series([2, 4, 6, 7, 8, 11, 13, 15, 18, 20])
}

print "...copying and renaming the columns to be PANAS_Q01...PANAS_Q21"
raw_cols = filter_list("PANAS", rawdf.columns)
panas_cols = new_cols_starting_with("PANAS", len(raw_cols))
df[panas_cols] = rawdf[raw_cols]

print "...cleaning up the item scores (ugh)"
def extract_number(score):
    res = re.search("[1-9]", score)
    if res: return int(res.group())
    else:   return 1    # note: 'not at all' gives nothing but should be 1
df[panas_cols] = df[panas_cols].applymap(extract_number)

print "...getting sub-scales"
sdf = compute_subscale_scores("PANAS", df[panas_cols], panas_scales)
df = append_subscales(df, sdf)

print "...getting total"
df["PANAS_Total"] = df[panas_cols].apply(sum, 1)


###
# AIM: Affect Intensity Measure
###

print "\nAIM"

aim_recoding = {"Never": 1, "Almost Never": 2, "Occasionally": 3, 
                "Usually": 4, "Almost Always": 5, "Always": 6}

# Questions to be negatively-keyed
aim_negatively_keyed = Series([12, 16, 18, 24, 26, 28, 29, 31, 33, 37, 40])

print "...copying and renaming the columns to be AIM_Q01...AIM_Q40"
raw_cols = filter_list("AIM", rawdf.columns)
aim_cols = new_cols_starting_with("AIM", len(raw_cols))
df[aim_cols] = rawdf[raw_cols]

print "...changing the text responses to a hard solid number"
df[aim_cols] = df[aim_cols].applymap(lambda x: aim_recoding[x])

print "...reversing negatively-keyed items"
subcols = aim_cols[aim_negatively_keyed-1].tolist()
df[subcols] = 7 - df[subcols]

print "...getting total"
df["AIM_Total"] = df[aim_cols].apply(sum, 1)


###
# Save!
###

print "\nSaving"

df.to_csv(path.join(behavdir, "ccd_all_filtered.csv"))
df.to_excel(path.join(behavdir, "ccd_all_filtered.xlsx"))

cols = ["study_id", "Age", "Sex"] + filter_list("Total", df.columns).tolist()
new_cols = [ x.replace("_Total", "") for x in cols ]
totals = df[cols]
totals.columns = new_cols
totals.to_csv(path.join(behavdir, "ccd_totals.csv"))
totals.to_excel(path.join(behavdir, "ccd_totals.xlsx"))

#widgets = ['Progress: ', Percentage(), ' ', Bar(), ' ', Counter(), '/%i' % len(details), ' ', ETA()]
#pb = ProgressBar(widgets=widgets, maxval=len(details)).start()
#pb.update(i)??
#pb.finish()??
