#!/usr/bin/env python

import random
import csv
csv_file_path = "./merged_oc_phen.tsv"
with open(csv_file_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t' )
    data = list(reader)
#column header 
data[0].append('Random')

#add random values 0 or 1 
for i in range(1, len(data)):
    data[i].append(random.choice([0, 1]))

#write into output
with open(csv_file_path, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerows(data)