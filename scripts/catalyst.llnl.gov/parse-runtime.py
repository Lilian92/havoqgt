#!/usr/bin/env python3

import glob
import csv

# Read in the raw outputs
results = []
for f in glob.glob('./output/output-pattern-matching-*.out'):
    print("Extracting data from: ", f)
    with open(f) as fi:
        cur = {}
        cur['LCC times'] = 0
        cur['LCC total runtime'] = 0
        cur['TPT times'] = 0
        cur['TPT total runtime'] = 0
        cur['TP times'] = 0
        cur['TP total runtime'] = 0
        for l in fi:
            if l.startswith("Running:"):
                # Extract the running parameters
                ls = l.split()
                cur['graph'] = ls[5]
                cur['scale'] = int(ls[5].split('rmat')[-1])
                cur['uniform random range'] = 0
                if int(ls[7]) == 0:
                    cur['degree_label'] = True
                else:
                    cur['degree based label'] = False
                    cur['uniform random range'] = int(ls[7])
                cur['nodes'] = 2
                cur['tasks per nodes'] = 4
                cur['pattern'] = ls[9]
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
 
        results.append(cur)


#write out to the CSV

with open('pattern-matching.csv', 'w') as f:
    w = csv.DictWriter(f, results[0].keys())
    w.writeheader()
    w.writerows(results)
