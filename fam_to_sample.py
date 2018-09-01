#!/usr/bin/env python

import os
import sys

"""Convert .fam files to .sample files for ADGC2.0 project.
.fam
ACT_ACT150544 ACT150544 0 0 -9 -9
ACT_ACT94662 ACT94662 0 0 1 1
ACT_ACT94442 ACT94442 0 0 1 1
ACT_ACT92870 38013785 0 0 -9 -9
ACT_ACT90327 ACT90327 0 0 1 1

to .sample
ID sex case
0 D B
sample1 M control
sample2 F control
sample3 F control
sample4 F case
sample4 NA NA
"""
USAGE = "./fam_to_sample.py fam_file"
HEADER = "ID sex case\n"
HEADER_DESC = "0 D B\n"
if len(sys.argv) != 2:
    print("USAGE: {}".format(USAGE))
    sys.exit(-1)

fam_file = sys.argv[1]
if not fam_file.endswith('.fam'):
    print("Input fam file must end with .fam!")
    sys.exit(-1)
sample_file = fam_file.replace('.fam', '.sample')
if os.path.exists(sample_file):
    print("{} already exists!".format(sample_file))
    sys.exit(-1)
print("Writing sample file: {}".format(sample_file))
with open(fam_file) as fam, open(sample_file, 'w+') as sample:
    # Write header lines
    sample.write(HEADER)
    sample.write(HEADER_DESC)
    for line in fam:
        sline = line.split()
        # 1. Combine column 1 and 2 with 3 dashes -> 1 first column
        sample_id = "{}__{}".format(sline[0], sline[1])
        # 2. Grab sex (5th column) converting 1 -> M and 2 -> F and 0 -> NA
        sex = sline[4]
        sample_sex = ""
        if sex == "1":
            sample_sex = "M"
        elif sex == "2":
            sample_sex = "F"
        elif sex == "-9":
            sample_sex = "NA"
        else:
            print("Error reading fam file. Unkown sex code: {}".format(sex))
            sys.exit(-1)
        # 3. Grab case/control status 6th column
        # Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
        pheno = sline[5]
        sample_pheno = ""
        if pheno == "1":
            sample_pheno = "control"
        elif pheno == "2":
            sample_pheno = "case"
        elif pheno == "-9":
            sample_pheno = "NA"
        else:
            print("Error reading fam file. Unknown phenotype code: {}".format(pheno))

        # Write data to new sample file
        sample.write("{} {} {}\n".format(sample_id, sample_sex, sample_pheno))
