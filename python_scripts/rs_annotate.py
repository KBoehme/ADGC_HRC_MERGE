outputfile = open('output.vcf', 'w+')

#read the refile into a dictionary
rsRefs = {}
prevPos = ''
thisPos = ''
chromNum = ''
with open('chrom21_ref.txt') as openRef:
	for line in openRef:
		words = line.split()
		chromNum = words[0]
		thisPos = words[1]
		if thisPos == prevPos:
			rsRefs[words[1]].update({(words[3] + words[4]):words[2]})
		else:
			rsRefs[words[1]] = {(words[3] + words[4]):words[2]}
		prevPos = thisPos

#open the vcf to be annotated and output into the open outputfile
with open('period_chr21.vcf') as openfile:
	for line in openfile:
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
		outputfile.write(line)
outputfile.close()
exit()
