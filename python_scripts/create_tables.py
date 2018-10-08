#!/usr/bin/env python
import fire
import pathlib

""" Script to generate tables for final documentation """

"""
Dataset & SNPs & Individuals & Size (GB) & Low Info ($<$.3) & Low Info (\%) \\ \hline
ACT1 & 39117111 & 2548 & 38 & 18569040 & 0.47 \\
ACT2 & 39127694 & 383 & 7.9 & 25246796 & 0.65 \\
ADC1 & 39117110 & 2738 & 43 & 11841628 & 0.30 \\
ADC2 & 39117110 & 927 & 16 & 22221054 & 0.57 \\
ADC3 & 39127707 & 1756 & 25 & 19351549 & 0.49 \\
ADC4 & 39127695 & 1048 & 17 & 19763571 & 0.51 \\
ADC5 & 39127695 & 1222 & 19 & 21052435 & 0.54 \\
ADC6 & 39127693 & 1326 & 20 & 20595016 & 0.53 \\
ADC7 & 39127693 & 1465 & 22 & 20351107 & 0.52 \\
ADNI & 39127689 & 691 & 13 & 23363694 & 0.60 \\
BIOCARD & 39127693 & 201 & 5.5 & 25338180 & 0.65 \\
CHAP2 & 39127694 & 743 & 13 & 22645052 & 0.58 \\
EAS & 39127696 & 283 & 6.3 & 25983088 & 0.66 \\
GSK & 39117177 & 1572 & 37 & 23043460 & 0.59 \\
LOAD & 39127741 & 4392 & 56 & 17215964 & 0.44 \\
MAYO & 39117108 & 1970 & 38 & 20315599 & 0.52 \\
MIRAGE & 39117109 & 1487 & 31 & 21069253 & 0.54 \\
NBB & 39127697 & 299 & 6.5 & 26732501 & 0.68 \\
OHSU & 39127684 & 607 & 14 & 25097262 & 0.64 \\
RMAYO & 39127694 & 428 & 8.3 & 25698016 & 0.66 \\
ROSMAP1 & 39065866 & 1649 & 29 & 21862619 & 0.56 \\
ROSMAP2 & 39127699 & 532 & 11 & 17278048 & 0.44 \\
TARC1 & 39127689 & 617 & 13 & 24333366 & 0.62 \\
TGEN2 & 39117117 & 1485 & 26 & 21011080 & 0.54 \\
UKS & 39117115 & 1741 & 28 & 20608444 & 0.53 \\
UMVUMSSM & 117372659 & 2467 & 36 & 66392871 & 0.57 \\
UMVUTARC2 & 39127693 & 540 & 11 & 21328491 & 0.55 \\
UPITT & 39127781 & 2212 & 32 & 18860900 & 0.48 \\
WASHU1 & 39127700 & 668 & 13 & 23744067 & 0.61 \\
WASHU2 & 37200188 & 235 & 6.4 & 22648158 & 0.61 \\
WHICAP & 39127693 & 644 & 11 & 23774868 & 0.61 \\ \hline
Totals & 1289129690 & 38876 & 652.9 & 717337177 & 0.56 \\ \hline
"""
from pathlib import Path
import subprocess
import os

def count_lines(filename):
    # return sum(1 for line in open(filename))
    return int(subprocess.check_output(["wc", "-l", filename], universal_newlines=True).split()[0])

def create_initial_data_summary(path_to_data, path_to_prune_snps):
    data = {}
    orig_data = Path(path_to_data)
    filter_data = Path(path_to_prune_snps)

    print("Working on variant counts")
    for fname in map(pathlib.Path, list( orig_data.glob('*/*info*.txt'))):
        study_name = fname.parts[-2]
        data[study_name] = {'var_count': count_lines(fname)}
    print("Working on sample counts")
    for fname in map(pathlib.Path, list( orig_data.glob('*/*.fam'))):
        study_name = fname.parts[-2]
        data[study_name]['sample_count'] = count_lines(fname)
    print("Working on low info counts")
    for fname in map(pathlib.Path, list( filter_data.glob('*/filters/*.txt'))):
        study_name = fname.parts[8]
        data[study_name]['low_info_count'] = count_lines(fname)
        data[study_name]['low_info_percent'] = data[study_name]['low_info_count'] / data[study_name]['var_count']

    print("""Dataset & SNPs & Individuals & Size (GB) & Low Info ($<$.3) & Low Info (\%) \\ \hline""")
    total_variants = 0
    total_sampls = 0
    total_low_info = 0
    info_pct_list = []
    for i,(k,v) in enumerate(data.items()):
        total_variants += v['var_count']
        total_sampls += v['sample_count']
        total_low_info += v['low_info_count']
        info_pct_list.append(v['low_info_percent'])

        vars = "{:,}".format(v['var_count'])
        sampls = "{:,}".format(v['sample_count'])
        lw_info = "{:,}".format(v['low_info_count'])
        info_prc = "{0:.2f}".format(v['low_info_percent'])
        out_string = "{} & {} & {} & {} & {} \\\\".format(k, vars, sampls, lw_info, info_prc)
        if i == len(data.keys()):
            out_string += " \\hline"
        print(out_string)
    total_variants = "{:,}".format(total_variants)
    total_sampls = "{:,}".format(total_sampls)
    total_low_info = "{:,}".format(total_low_info)
    avg_low_info = sum(info_pct_list) / float(len(info_pct_list))
    avg_low_info = "{0:.2f}".format(avg_low_info)

    print("Totals & {} & {} & {} & {} \\\\ \\hline".format(total_variants, total_sampls, total_low_info, avg_low_info))
    print(data)


def create_final_sum_stats():
    """ADGC HRC (Full) 	&	ADGC HRC (Unrelated) \\ \hline
    row 1 = snps
    row 2 = samples
    row 3 = genotyping rate
    """
    pass

def create_final_full_sum():
    """Study 	&	Sex(M/F)	&	Cases/Controls	&	Sample Size	\\ \hline"""
    pass

def create_final_unrelated_sum():
    """Study 	&	Sex(M/F)	&	Cases/Controls	&	Sample Size	\\ \hline
    ACT 1   1,090/1,458 567/1,701/280   2,548
    """
    pass

if __name__ == '__main__':
  fire.Fire()
