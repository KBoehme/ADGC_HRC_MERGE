#!/usr/bin/env python

""" Script to compare the headers of different covariant files. """

import sys
# import pandas as pd
import fire
from pprint import pprint
from pathlib import Path

def combine_covars(covar_file_dir, combined_output):
    cvd = Path(covar_file_dir)
    covar_files = cvd.glob("*.txt")
    my_headers_dict = {}
    my_headers = []
    header_to_use = []
    from collections import defaultdict
    total_samples = defaultdict(int)
    for covar in covar_files:
        print(type(covar))
        with open(covar, 'r') as c:
            header = next(c)
            header_list = header.strip().split()
            header_set = set(header_list)
            print(f"file={covar},header size={len(header_list)} ")
            set(header.strip().split())
            if covar.name == "umvumssm_a.covar.2013Dec06.txt":
                header_to_use = header_list
            print(f"file={covar},hlist={len(header_list)},hset={len(header_set)}")
            my_headers_dict[covar] = header_set
            my_headers.append(header_set)

    u = set.intersection(*my_headers)
    header_union = set.union(*my_headers)
    print(f"intersection={u}")
    for covar, header in my_headers_dict.items():
        diff = list(header-u)
        print(f"covar={covar},unique_to_covar={diff}")

    print(f"header_union={header_union}")
    print(f"header_to_use={header_to_use}")
    diff = set(header_to_use) - header_union
    print(f"Difference between header to use and header union: {diff}")
    if len(diff) > 0:
        print("Problem, header to use doesnt contain the entire union of fields.")
        sys.exit(1)

    MISSING = "-9"
    with open(combined_output, 'w') as outfile:
        outfile.write('\t'.join(header_to_use)+"\n")
        for covar, header in my_headers_dict.items():
            print(f"On covar={covar}")
            with open(covar) as c:
                my_header = ""
                for i, line in enumerate(c):
                    if i == 0:
                        my_header = line.strip().split()
                        print(my_header)
                    else:
                        total_samples[covar] += 1
                        outline = []
                        sline = line.strip().split()
                        zip_line = dict(zip(my_header, sline))
                        for item in header_to_use:
                            outline.append(zip_line.get(item, MISSING))
                        outfile.write('\t'.join(outline)+"\n")
    print("Total samples found: {}".format(sum(total_samples.values())))
    print(pprint(total_samples))

if __name__ == '__main__':
  fire.Fire()
