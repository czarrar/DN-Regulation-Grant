#!/usr/bin/env bash

njobs=15

cd /home/data/Projects/CCD/scripts/01_preprocessing

# Preprocess Functional
cat ../../msit_subject_list.txt | parallel -j $njobs --eta './01_preprocess_func.rb -o MSIT -n msit -s {}'

# Register Functional
cat ../../msit_subject_list.txt | parallel -j $njobs --eta './03_register_func.rb -n msit -r 1 2 -s {}'

# Regress Nuisance
cat ../../msit_subject_list.txt | parallel -j $njobs --eta './04_nuisance_func.rb -n msit -r 1 2 -s {}'

# Apply Registration
cat ../../msit_subject_list.txt | parallel -j $njobs --eta './05_applyreg_func.rb -n msit -r 1 2 -s {}'
