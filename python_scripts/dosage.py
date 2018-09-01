#!/usr/bin/env python
""" Convert vcf to dosage. Requires AF in order to differentiate major/minor reference allele. If the AF is > .5 then
the reference is taken as the minor allele. If the AF is lower than .5 the ref is taken as the major allele.
"""
import sys
import re

in_vcf = sys.argv[1]
out_vcf = sys.argv[2]
CHUNK_SIZE = 1000

def write_chunk_to_file(out_file, data):
	""" Write data to outfile. data list list and will be written to file
	line by line
	"""
	for entry in data:
		out_file.write(entry+"\n")

with open(in_vcf) as vcf, open(out_vcf, 'w') as out_vcf:
	##FORMAT=<ID=GP,Type=Float,Number=G,Description="Genotype call probabilities">
	##FORMAT=<ID=GT,Type=String,Number=1,Description="Genotype call probabilities, threshholded at 0.51">
	##FORMAT=<ID=DS,Type=String,Description="Genotype dosage.">
	##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
	#CHROM   POS             ID      		 REF     ALT     QUAL    FILTER  INFO    FORMAT  sample_0
	#22      16050435        22:16050435,22  T       C       .       .       AF=1       GT:GP   0/0:0.999,0.001,0
	#chr3 186552857 . A G . . . GT:DP:DS 0/1:0,0.995,0.005:1.005 0/0:0.8,0.2,0:0.2
	line_list = []
	chunk_number = 0
	for i, line in enumerate(vcf.readlines()):
		sline = line.split()
		if line.startswith("##"):
			if line.startswith("##INFO"): # Stick in new format field.
				out_vcf.write("##FORMAT=<ID=DS,Type=String,Description=\"Genotype dosage.\">\n")
			out_vcf.write(line)
			pass
		elif line.startswith("#"):
			sline = line.lstrip("#").split()
			out_vcf.write(line)
			header = sline
		else:
			hline = dict(zip(header, sline))
			af = hline['INFO'].split('=')
			m = re.search("AF=(\d+)", hline['INFO'])
			af = int(m.group(1))
			lines_to_write = []

			# Use current line, just add :DS to end of FORMAT column.
			current_line = f"{hline['CHROM']}\t{hline['ID']}\t{hline['REF']}\t{hline['ALT']}\t{hline['QUAL']}\t{hline['FILTER']}\t{hline['INFO']}\t{hline['FORMAT']}:DS"
			all_sample_data = []
			for sample in sline[9:]:
				format_list = hline['FORMAT'].split(":")
				hsample_data = dict(zip(format_list, sample.split(":")))
				gps = hsample_data['GP'].split(',')
				gps = list(map(float, gps))
				# print(gps)
				if len(gps) != 3:
					print("Expected 3 GP values found: {}".format(len(gps)))
					sys.exit(1)
				ds = 0
				# print(gps)
				if af >= 0.5:
				 	ds = gps[0] + gps[0] + gps[1]
				else:
					ds = gps[2] + gps[2] + gps[1]
				# Use current sample data, just add a :<calculated dosage> to the end.
				all_sample_data.append(f"{sample}:{ds}")
			joined_sample_data = '\t'.join(all_sample_data)
			current_line = f"{current_line}\t{joined_sample_data}"
			line_list.append(current_line)
			len_line_list = len(line_list)
			if len_line_list > 1 and len_line_list % CHUNK_SIZE == 0:
				# Write chunk to file
				chunk_number += 1
				print("Writing samples: {}".format(CHUNK_SIZE*chunk_number))
				write_chunk_to_file(out_vcf, line_list)
				all_sample_data = []
	if line_list:
		write_chunk_to_file(out_vcf, line_list)
