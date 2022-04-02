#!/usr/bin/python3.8
 # -*- coding: utf-8 -*-

import csv
import subprocess
import sys


region_file=sys.argv[1]
file_handle=sys.argv[2]
outfile="./Output/coverage_dataframe.tsv"

with open(region_file) as l:
	r=l.readlines()
region_list=[el.strip() for el in r]

with open(file_handle) as l:
	h=l.readlines()
file_list=[el.strip() for el in h]

with open(outfile, 'w') as write_obj:
	csv_writer = csv.writer(write_obj, dialect="excel", delimiter='\t')
	header=["Sample_ID","Region","Length (bp)","E[DoC]","BoC"]
	csv_writer.writerow(header)
	for file_in in file_list:
		region_index=0
		for region in region_list:
			region_index+=1
			region_values=subprocess.run(["grep",region,file_in], text=True, stdout=subprocess.PIPE)
			csv_lines=region_values.stdout.split("\n")
			csv_reader = csv.reader(csv_lines[:-1], dialect="excel", delimiter='\t')
			mean_DoC=0.0
			BoC=1.0
			leg=0
			Doc=0.0
			for row in csv_reader:
				mean_DoC+=float(row[3])*float(row[4])/float(row[5])
				leg=int(row[5])
				if float(row[3])==0:
					BoC=round((1-float(row[6])),3)
			DoC=round(mean_DoC,3)
			outline=[file_in,"R"+str(region_index),leg,DoC,BoC]
			csv_writer.writerow(outline)
