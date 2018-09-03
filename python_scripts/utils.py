""" Collection of utility scripts"""

def get_size_file_mb(filename):
    return int(os.path.getsize(filename)/(1024*1024))

def get_mem(wildcards):
    print(wildcards)
    return int(get_size_mb('output.dat')*10)

def rs_annotate(input_vcf, rs_reference, chrom, output):
    import gzip

    outputfile = gzip.open(output, 'w')

    #read the refile into a dictionary
    rsRefs = {}
    prevPos = ''
    thisPos = ''
    chromNum = ''
    with gzip.open(rs_reference, 'r') as openRef:
        for line in openRef:
            line = line.decode('utf-8')
            words = line.split()
            chromNum = words[0]
            thisPos = words[1]
            if thisPos == prevPos:
                rsRefs[words[1]].update({(words[3] + words[4]):words[2]})
            else:
                rsRefs[words[1]] = {(words[3] + words[4]):words[2]}
            prevPos = thisPos

    #open the vcf to be annotated and output into the open outputfile
    CHUNK = 10000
    with gzip.open(input_vcf, 'r') as openfile:
        for i, line in enumerate(openfile):
            line = line.decode('utf-8')
            if i % CHUNK == 0:
                print("Working on line number: {}".format(i))
            wrds = line.split()
            if wrds[0][0] != '#':
                val = rsRefs.get(wrds[1], None)
                position = line.find(wrds[2])
                if val is not None:
                    val2 = val.get((wrds[3] + wrds[4]), None)
                    if val2 is not None:
                        line = line[0:position] + val2 + line[(position + len(wrds[2])):len(line)]
                    else:
                        line = line[0:position] + wrds[2][0:wrds[2].find(',')] + line[(position + len(wrds[2])):len(line)]
                else:
                    line = line[0:position] + wrds[2][0:wrds[2].find(',')] + line[(position + len(wrds[2])):len(line)]
                line = chromNum + line[(line.find(wrds[0]) + len(wrds[0])):len(line)]
            outputfile.write(line.encode())
            if i == 0:
                outputfile.write("##contig=<ID={}>\n".format(chrom).encode())
    outputfile.close()


def sample_from_fam(input, output):
    print("Running create_sample_file_from_fam with input={} and output={}".format(input, output))
    import os
    import sys
    """Convert .fam files to .sample files for ADGC2.0 project.
    .fam

    Family ID ('FID')
    Within-family ID ('IID'; cannot be '0')
    Within-family ID of father ('0' if father isn't in dataset)
    Within-family ID of mother ('0' if mother isn't in dataset)
    Sex code ('1' = male, '2' = female, '0' = unknown)
    Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

    ACT_ACT150544 ACT150544 0 0 -9 -9
    ACT_ACT94662 ACT94662 0 0 1 1
    ACT_ACT94442 ACT94442 0 0 1 1
    ACT_ACT92870 38013785 0 0 -9 -9
    ACT_ACT90327 ACT90327 0 0 1 1

    converted from 1/0 to 2/1 case/control coding).
    to .sample
    ID_1 ID_2 missing case cov_1
    0 0 0 B C
    sample1 1 0
    sample2 2 0
    sample3 2 0
    sample4 2 1
    sample4 0 NA
    """
    HEADER = "ID_1 ID_2 missing sex case\n"
    HEADER_DESC = "0 0 0 D B\n"
    if len(input) != 1 or len(output) != 1:
        sys.exit(-1)

    fam_file = input[0]
    sample_file = output[0]
    print("Writing to output file: {}".format(sample_file))
    if not fam_file.endswith('.fam'):
        print("Input fam file must end with .fam!")
        sys.exit(-1)
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

            sample_unique_id = "{}___{}".format(sline[0], sline[1])
            sex = sline[4]
            sample_sex = ""
            if sex == "1":
                sample_sex = "1"
            elif sex == "2":
                sample_sex = "2"
            elif sex == "-9":
                sample_sex = "0"
            else:
                print("Error reading fam file. Unkown sex code: {}".format(sex))
                sys.exit(-1)
            # 3. Grab case/control status 6th column
            # Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            pheno = sline[5]
            sample_pheno = ""
            if pheno == "1":
                sample_pheno = "0"
            elif pheno == "2":
                sample_pheno = "1"
            elif pheno == "-9":
                sample_pheno = "NA"
            else:
                print("Error reading fam file. Unknown phenotype code: {}".format(pheno))
                sys.exit(-1)

            # Write data to new sample file
            # Use a combined ID_1___ID_2 as the first ID because thats sane.
            sample.write("{} {} 0 {} {}\n".format(sample_unique_id, sline[1], sample_sex, sample_pheno))
