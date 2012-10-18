#!/usr/bin/env bash

njobs=15

cd /home/data/Projects/CCD/scripts/01_preprocessing

# Preprocess Functional
#cat ../../subject_list.txt | parallel -j $njobs --eta './01_preprocess_func.rb -o TRAIN -n rest -s {}'

# Register Functional
cat ../../subject_list.txt | parallel -j $njobs --eta './03_register_func.rb -n rest -r 1 -s {}'

# Regress Nuisance
cat ../../subject_list.txt | parallel -j $njobs --eta './04_nuisance_func.rb -n rest -r 1 -s {}'

# Apply Registration
cat ../../subject_list.txt | parallel -j $njobs --eta './05_applyreg_func.rb -n rest -r 1 -s {}'
