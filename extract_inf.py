#/bin/python
# -*- coding:utf-8 -*-
# to extract overexpressed reads from fastqc results
# demo :  python extract_inf.py  -i fastqc_data.txt -o output_file

import sys, os, re, optparse

optParser = optparse.OptionParser(
        usage = "python %prog [options]",
        description = "This script calls .",
        epilog = "Example: python extract_inf.py  -i fastqc_data.txt -o output_file" )

optParser.add_option( "-i", "--in", type="string", dest="input", help = "input fasta file [REQUIRED]" )
optParser.add_option( "-o", "--out", type="string", dest="output", help = "output fasta file [REQUIRED]" )

(opts, args) = optParser.parse_args()

if not (opts.input or opts.output):
        optParser.print_help()
        sys.exit(1)


file_r = open(opts.input,"r")
file_w = open(opts.output,"w")
OverexpressRecords = []
for line in file_r:
	if line[0:4] == ">>Ov":
		for line in file_r:
			if line[0] != ">":
				OverexpressRecords.append(line)
			else:
				break
		file_w.writelines(OverexpressRecords)
		break
file_r.close()
file_w.close()

