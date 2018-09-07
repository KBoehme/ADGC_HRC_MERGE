#!/usr/bin/env python

import sys
from collections import defaultdict

"""
1       324822  A       T       0.0733954
1       324822  A       C       6.49421e-05
1       324822  A       G       2.16474e-05


1       324822  A       T       0.0733954
1       324822  A       C       6.49421e-05
1       324822  A       G       2.16474e-05
"""

def write_low_info_and_dups_filter_file(info_file, gnomad_af_file, low_info_threshold, output, dups_output):
    """ Python program that takes input list of snps and filters on both low
    info and duplicate variants that are more rare given gnomad dataset.

    If remove_dups=False, this will simply write rarer multi-alts to a file."""

    print(info_file)
    print(gnomad_af_file)
    print(low_info_threshold)
    print(output)
    # Write header to output file.
    # Important that its space delimited and the first 4 columns are there.
    outfile = open(output, 'w+')
    dups_outfile = open(dups_output, 'w+')

    outfile.write("SNPID rsid chromosome position ref alt\n")

    # Open input and loop through looking for duplicates
    in_snps = defaultdict(list)
    in_snps_maf_dict = {}
    with open(info_file, 'r') as f:
        #SNP REF(0) ALT(1) ALT_Frq MAF AvgCall Rsq
        # 1:13380 C G 0.00003 0.00003 0.99997 0.00022
        next(f) # Skip header line, Nice because no check every iteration
        for line in f:
            # Loop through file, writing low info variants
            sline = line.split()
            chrom_pos = sline[0]
            chr, pos = chrom_pos.split(":")
            ref = sline[1]
            alt = sline[2]
            maf = float(sline[4])
            # Found one with -
            if sline[6] == "-" or float(sline[6]) < float(low_info_threshold):
                out_line = [chr, chrom_pos, chr, pos, ref, alt]
                # print("Writing low info variant: {}".format(' '.join(out_line)))
                outfile.write(" ".join(out_line)+"\n")
            else:
                key = "{}+{}+{}".format(chr, pos, ref)
                in_snps_maf_dict["{}+{}".format(key, alt)] = maf
                in_snps[key].append(alt)

    # Read in the gnomad data
    # gnomad_af = {}
    # with open(gnomad_af_file, 'r') as f:
    #     for line in f:
    #         sline = line.split()
    #         if len(sline[2]) > 1 or len(sline[3]) > 1: # Skip anything not a snp
    #             pass
    #         else:
    #             #1       12272   G       A       0.0
    #             key = "{}+{}+{}+{}".format(sline[0], sline[1], sline[2], sline[3])
    #             gnomad_af[key] = float(sline[4])

    # Now deal with duplicates
    missings = 0
    found = 0
    for k,v in in_snps.items():
        chr, pos, ref = k.split('+')
        if len(v) > 1: # duplicate
            # Resolve duplicates
            # print("Duplicate variant found: key= {}, alts = {}".format(k, v))
            highest_maf = -1
            highest_maf_alt = None
            for alt_dup in v:
                variant_key = "{}+{}".format(k, alt_dup)
                missings += 1
                af = in_snps_maf_dict[variant_key]
                if af is None:
                    raise RuntimeError("Didnt find a MAF for variant: {}".format(variant_key))
                if af >= highest_maf:
                    if highest_maf_alt is not None: # Write current one to remove file
                        #SNPID rsid chromosome position ref alt
                        # 1:13380 1:13380 1 13380 C G
                        # 1:54676 1:54676 1 54676 C T
                        chrom_pos = "{}:{}".format(chr, pos)
                        out_line = [chrom_pos,
                            chrom_pos,
                            chr,
                            pos,
                            ref,
                            highest_maf_alt] # Must write previously highest maf alt
                        # print("Writing variant: {}".format(' '.join(out_line)))
                        dups_outfile.write(" ".join(out_line)+"\n")
                    highest_maf_alt = alt_dup # Now update highest maf allele
                    highest_maf = af
                else:
                    if highest_maf is not None:
                        chrom_pos = "{}:{}".format(chr, pos)
                        out_line = [chrom_pos,
                            chrom_pos,
                            chr,
                            pos,
                            ref,
                            alt_dup] # Must write current alt
                        # print("Writing variant: {}".format(' '.join(out_line)))
                        outfile.write(" ".join(out_line)+"\n")
                    # It wasnt as large as highest_maf.
    print("Done writing {}".format(output))
    print("Missings = {}".format(missings))
    print("Found in gnomad = {}".format(found))
    outfile.close()
    dups_outfile.close()

write_low_info_and_dups_filter_file(*sys.argv[1:])
