#!/usr/bin/env python

from cliutils import Process
from multiprocessing import Pool
from glob import glob
from os import path
from pandas import *
import numpy as np

basepath = "/home2/data/Projects/CCD"

###
# CCB: Generate paths to preprocessed MSIT functionals
###

preproc = path.join(basepath, "CCB/preprocessed")
ccb_filt = read_csv(path.join(basepath, "scripts/data/ccb_subject_info_filt.csv"))
ccb_filt = ccb_filt[ccb_filt.condition=="MSIT"]
ccb_filt = ccb_filt[["subject", "scan", "run"]]

def gen_inpaths(row):
    find_path = path.join(preproc, 
        "%s_scan%i_*" % (row['subject'], row['scan']), 
        "swnmrda%s_scan%i_MSIT%i_*_forfsl.nii.gz" % (row['subject'], row['scan'], row['run'])
    )
    paths = glob(find_path)
    if len(paths) == 0:
        print "Error: no paths found for %s" % find_path
        print row
        return (None, None)
    if len(paths) > 1:
        print "Error: %i paths found for %s" % (len(paths), find_path)
        print row
    inpath = paths[0].rstrip(".nii.gz")
    return inpath

def gen_outpaths(row):
    outpath = path.join(basepath, "analysis/subjects", row['subject'], 
                            "msit_%02i" % row['scan'], 
                            "run_%02i" % row['run'], 
                            "task_hrf.feat")
    return outpath

ccb_inpaths = ccb_filt.apply(gen_inpaths, axis=1)
ccb_outpaths = ccb_filt.apply(gen_outpaths, axis=1)

args = zip(ccb_inpaths.tolist(), ccb_outpaths.tolist())
    
def run_process(arg):
    from cliutils import Process
    from os import path
    
    template = "level1_ccb_template.fsf"
    template_input = "/home/data/Projects/CCD/CCB/preprocessed/CCB002_scan1_20100602/swnmrdaCCB002_scan1_MSIT1_20100602_forfsl"
    template_output = "/home2/data/Projects/CCD/analysis/subjects/CCB002/msit1/run_01/task_hrf.feat"
    
    row = {
        'input': arg[0], 
        'output': arg[1]
    }
    
    basefsf = row['output'].rstrip('.feat') + 'f'
    outfsf = basefsf + '.fsf'
    
    cmd = "mkdir -p %s" % path.dirname(basefsf)
    print cmd
    print Process(cmd)
    
    print '%s => %s' % (template, outfsf)
    f = open(template, 'r')
    text = f.read()
    text = text.replace(template_input, row['input'])
    text = text.replace(template_output, row['output'])
    f = open(outfsf, 'w')
    f.write(text)
    f.close()
    
    cmd = "feat_model %s" % basefsf
    print cmd
    print Process(cmd)
    
    cmd = "feat %s" % outfsf
    print cmd
    print Process(cmd)
    
    regdir = path.join(row['output'], 'reg')
    cmd = "mkdir %s" % regdir
    print cmd
    print Process(cmd)
    
    cmd = "cp /home2/data/PublicProgram/fsl-4.1.9/etc/flirtsch/ident.mat %s" % path.join(regdir, "example_func2standard.mat")
    print cmd
    print Process(cmd)
    
    cmd = "cp /home2/data/PublicProgram/fsl-4.1.9/etc/flirtsch/ident.mat %s" % path.join(regdir, "highres2standard.mat")
    print cmd
    print Process(cmd)
    
    cmd = "ln -sf /home2/data/PublicProgram/fsl-4.1.9/data/standard/MNI152_T1_4mm_brain.nii.gz %s" % path.join(regdir, "highres.nii.gz")
    print cmd
    print Process(cmd)

    cmd = "ln -sf /home2/data/PublicProgram/fsl-4.1.9/data/standard/MNI152_T1_4mm_brain.nii.gz %s" % path.join(regdir, "standard.nii.gz")
    print cmd
    print Process(cmd)

p = Pool(35)
p.map(run_process, args)
