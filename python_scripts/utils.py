""" Collection of utility scripts"""
import gzip

def get_size_file_mb(filename):
    return int(os.path.getsize(filename)/(1024*1024))

def get_mem(wildcards):
    print(wildcards)
    return int(get_size_mb('output.dat')*10)


def read_ref_file(rs_reference):
    rsRefs = {}
    with gzip.open(rs_reference, 'rt') as openRef:
        for line in openRef:
            line = line
            words = line.split()
            chrom = words[0]
            pos = words[1]
            rsid = words[2]
            ref = words[3]
            alt = words[4]
            rsRefs[f"{chrom}+{pos}+{ref}+{alt}"] = rsid
    print("Done reading ref file")
    return rsRefs

def create_update_map(gen_file, ref_file, update_map):
    """ Creates 12 column variant id map updater for qctools"""
    gen_file = gen_file
    ref_file = ref_file
    out_file = update_map

    rs_ref = read_ref_file(ref_file)
    with gzip.open(gen_file, 'rt') as openfile, open(out_file, 'w') as out:
        out.write("SNPID rsid chromosome position a1 a2 new_SNPID new_rsid new_chrom new_pos new_a1 new_a2\n")
        for line in openfile:
            sline = line.split()
            SNPID, rsid, chromosome, position, a1, a2 = sline[:6]
            outline = []
            outline.extend(sline[:6])
            #22 22 22:16050435 16050435 T C 0.999
            #22 22 rs544901529 16352459 C T
            key = f"{SNPID}+{position}+{a1}+{a2}"
            new_id = rs_ref.get(key, None)
            if new_id is not None:
                outline.extend([SNPID, new_id, chromosome.split(":")[0], position, a1, a2])
                if len(outline) != 12:
                    sys.exit("Failing outline not equal to 12.")
                out.write(' '.join(outline)+"\n")

def rs_update_gen(gen_file, ref_file, output):
    """ Loop through gen_file and update rsid column to real rsid.

    SNP ID, RS ID of the SNP, base-pair position of the SNP, the allele coded A and the allele coded B
    22      22                  22:17060707 17060707 G A

     """
    import gzip

    outputfile = gzip.open(output, 'wt')
    #read the refile into a dictionary
    rsRefs = read_ref_file(ref_file)

    #open the vcf to be annotated and output into the open outputfile
    CHUNK = 1000
    print("Start parsing input gen")
    with gzip.open(gen_file, 'rt') as openfile:
        write_chunk = []
        for i, line in enumerate(openfile):
            line = line
            if i != 0 and i % CHUNK == 0:
                print("Working on line number: {}".format(i))
                if write_chunk is not None:
                    outputfile.writelines(write_chunk)
                    write_chunk = []
            sline = line.split()
            chrom = sline[0]
            rsid_chrom = sline[1]
            chrom_pos = sline[2]
            pos = sline[3]
            ref = sline[4]
            alt = sline[5]
            key = f"{chrom}+{pos}+{ref}+{alt}"
            # If we cant find an rsid, create a detailed variant id
            detailed_variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            new_id = rsRefs.get(key, detailed_variant_id)
            outline = [chrom, new_id, pos, ref, alt]
            outline.extend(sline[6:]) # Add the rest of the line unchanged
            write_chunk.append(" ".join(outline)+"\n")
        if write_chunk:
            outputfile.writelines(write_chunk)
    outputfile.close()
    print("Done writing output gen: {}".format(output))


def rs_annotate(input_vcf, rs_reference, chromNum, output):
    import gzip

    outputfile = gzip.open(output, 'wt')
    chromNum = str(chromNum)
    #read the refile into a dictionary
    rsRefs = read_ref_file(ref_file)

    #open the vcf to be annotated and output into the open outputfile
    CHUNK = 1000
    print("start parsing input vcf")
    with gzip.open(input_vcf, 'rt') as openfile:
        write_chunk = []
        for i, line in enumerate(openfile):
            line = line
            if i != 0 and i % CHUNK == 0:
                print("Working on line number: {}".format(i))
                if write_chunk is not None:
                    outputfile.writelines(write_chunk)
                    write_chunk = []
            if not line.startswith('#'):
                sline = line.split()
                chrom = sline[0]
                pos = sline[1]
                id = sline[2]
                ref = sline[3]
                alt = sline[4]
                key = f"{chrom}+{pos}+{ref}+{alt}"
                new_id = rsRefs.get(key, id.replace(f",{chromNum}", ""))
                outline = [chromNum, pos]
                outline.append(new_id)
                outline.extend(sline[3:]) # Add the rest of the line unchanged
                outline += "\n"
                write_chunk.append("\t".join(outline))
            else:
                write_chunk.append(line)
                if i == 0:
                    write_chunk.append(f"##contig=<ID={chromNum}>\n")
        if write_chunk:
            outputfile.writelines(write_chunk)
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
