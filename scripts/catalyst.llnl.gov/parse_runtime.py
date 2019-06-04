#!/usr/bin/env python3

import glob
import csv
import math

# Read maximum degrees
max_degs = {}
for f in glob.glob('./output/output_gen_rmat_*.out'):
    print("Extracting data from: ", f)
    with open(f) as fi:
        scale = 0
        max_deg = 0
        for l in fi:
            if l.startswith("Running"):
                # Extract the running parameters
                ls = l.split()
                scale = int(ls[3])
            if l.startswith("Max Degree = "):
                # Extract running time of vertex data generating
                ls = l.split()
                max_deg = int(ls[3])
        max_degs[scale] = max_deg

# Read in the raw outputs
results = []
for f in glob.glob('./output/output_pattern_matching_*.out'):
    print("Extracting data from: ", f)
    with open(f) as fi:
        cur = {}
        cur['LCC times'] = 0
        cur['LCC total runtime'] = 0
        cur['TPT times'] = 0
        cur['TPT total runtime'] = 0
        cur['TP times'] = 0
        cur['TP total runtime'] = 0 
        cur['subpattern count'] = 'unknown'
        cur['Pattern Matching Time runtime'] = 'unknown'
        cur['max degree'] = 'unknown'
        cur['log2 of max degree'] = 'unknown'
        for l in fi:
            if l.startswith("Running:"):
                # Extract the running parameters
                ls = l.split()
                cur['graph'] = ls[5].split('/')[-1]
                cur['scale'] = int(ls[5].split('rmat')[-1])
                if cur['scale'] in max_degs:
                    cur['max degree'] = max_degs[cur['scale']]
                    cur['log2 of max degree'] = math.log2(cur['max degree'])
                else:
                    print('rmat', cur['scale'], ' are not found in generator')
                cur['uniform random range'] = 0
                if int(ls[7]) == 0:
                    cur['degree label'] = True
                else:
                    cur['degree label'] = False
                    cur['uniform random range'] = int(ls[7])
                cur['nodes'] = 2
                cur['tasks per nodes'] = 4
                cur['pattern'] = int(ls[9].split('/')[-1])
            if l.startswith("Pattern Matching Time | Vertex Data DB :"):
                # Extract running time of vertex data generating
                ls = l.split()
                cur['mate data generating time'] = float(ls[8])
            if l.startswith("Pattern Matching Time | Local Constraint Checking :"):
                # Extract running time of local constraint checking
                ls = l.split()
                cur['LCC times'] += 1
                cur['LCC total runtime'] += float(ls[8])
            if l.startswith("Pattern Matching Time | Token Passing (Traversal)"):
                # Extract running time of token passing traversal
                ls = l.split()
                cur['TPT times'] += 1
                cur['TPT total runtime'] += float(ls[9])
            if l.startswith("Pattern Matching Time | Token Passing ["):
                # Extract running time of token passing
                ls = l.split()
                cur['TP times'] += 1
                cur['TP total runtime'] += float(ls[8])
            if l.startswith("Pattern Matching Time | Pattern"):
                # Extract running time of pattern matching
                ls = l.split()
                cur['Pattern Matching Time runtime'] = float(ls[7])
            if l.find("Global Subpattern Count") != -1:
                # Extract number of matched subpattern
                ls = l.split()
                cur['subpattern count'] = int(ls[8])
 
        results.append(cur)


#write out to the CSV

with open('pattern_matching.csv', 'w') as f:
    w = csv.DictWriter(f, results[0].keys())
    w.writeheader()
    w.writerows(results)
