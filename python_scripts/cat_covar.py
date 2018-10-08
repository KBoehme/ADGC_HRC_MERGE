#!/usr/bin/env python

""" Script to compare the headers of different covariant files. """

import sys
# import pandas as pd
import fire
from pprint import pprint
from pathlib import Path


def read_idmap(shortid_map):
	original_to_small = {}
	small_to_original = {}
	with open(shortid_map) as idmapfile:
		for line in idmapfile:
			sline = line.strip().split()
			fid_id_key = "{}:{}".format(sline[2],sline[3])
			original_key =  "{}:{}".format(sline[0],sline[1])
			original_to_small[original_key] = fid_id_key
			small_to_original[fid_id_key] = original_key
	return original_to_small, small_to_original

def add_pcs(covar_file, pcs, shortid_map, output):
	"""
	#eigvals:    80.126    17.338    14.071    13.632    12.164     8.620     8.295     8.103     7.500     7.173
	F1:I1     0.0026     -0.0059      0.0056     -0.0049     -0.0131      0.0138      0.0024     -0.0023     -0.0114      0.0048             Case
	"""
	orig_id_map, short_id_map = read_idmap(shortid_map)
	pc_data = {}
	with open(pcs) as pcfile:
		for i,line in enumerate(pcfile):
			line = line.strip()
			sline = line.split()
			if line.startswith("#"):
				header = sline
			else:
				orig_key = short_id_map[sline[0]]
				pc_data[orig_key] = sline[1:11]
	print("Number of samples with PC data in dict = {}".format(len(pc_data.keys())))
	print(f"Number of lines in PC data file = {i}")
	print(type(pc_data))
	found_data = 0
	# Finally lets start writing the output file
	with open(covar_file) as covar, open(output, 'w') as out:
		header = next(covar).split()
		header.extend(["pc{}_merged".format(num) for num in range(1,11)])
		header.append("older_65_control")
		# print(header)
		out.write(' '.join(header)+'\n')
		for line in covar:
			sline = line.split()
			covar_fid, covar_iid = sline[0], sline[1]
			my_pcs = pc_data.get("{}:{}".format(covar_fid, covar_iid))
			pcs_to_add = []
			if my_pcs is not None:
				# We found something
				pcs_to_add = my_pcs
				found_data += 1
			else:
				pcs_to_add = ['-9' for i in range(0,10)]
			sline.extend(pcs_to_add)
			if sline[8] == "1" and (int(float(sline[13]) ) > 65 and sline[13] != "-9"):
				old_control = "1"
			else:
				old_control = "0"
			sline.append(old_control)
			out.write(' '.join(sline)+'\n')
	print("We matched pc data for {} samples".format(found_data))
	if found_data != len(pc_data.keys()):
		print("Failed to match up every sample with pc data!!!!!!!!!!!!!!!")
		sys.exit(1)

def get_umv_map(orig_id_map):
	umv_map = {}
	for k,v in orig_id_map.items():
		if k.startswith('0___'):
			fixed_k = k.replace('0___', 'UMVUMSSM_')
			split_id = fixed_k.split('_')
			fid = "{}_{}".format(split_id[0], split_id[1])
			iid = k.split("-")[-1]
			new_corrected_key = "{}___{}:{}".format(fid, iid, iid)
			umv_map[new_corrected_key] = k
	return umv_map

def update_ids(covar_file, idmap, output):
	orig_id_map, _ = read_idmap(idmap)
	umv_ids = get_umv_map(orig_id_map)
	found_data = 0
	# Finally lets start writing the output file
	with open(covar_file) as covar, open(output, 'w') as out:
		header = next(covar).split()
		out.write(' '.join(header)+'\n')
		for line in covar:
			sline = line.split()
			covar_fid, covar_iid = sline[0], sline[1]
			covar_key = "{}___{}:{}".format(covar_fid, covar_iid, covar_iid)
			if covar_key in orig_id_map:
				# Easy match.
				found_data += 1
				sline[0] = "{}___{}".format(covar_fid, covar_iid)
				sline[1] = covar_iid
			else:
				correct_id = umv_ids.get(covar_key)
				if correct_id is not None:
					found_data += 1
					sline[0], sline[1] = correct_id.split(":")
				else:
					print("FAILED TO FIND A MATCH!!!!!!!!!!!")
			out.write(' '.join(sline)+'\n')
		# Add final missing sample.
		missing_entry = ["ADC4_NACC101921___NACC101921", "NACC101921", "ADC4"]
		for i in range(len(missing_entry), len(header)):
			missing_entry.append("-9")
		out.write(' '.join(missing_entry)+'\n')
	print("We renamed ids for {} samples".format(found_data))
	if found_data != len(orig_id_map):
		print("Failed to rename all samples")


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
