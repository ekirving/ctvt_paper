#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, re

if len(sys.argv) != 3:
    print "Error: Missing FAM file"
    quit()

fam_path = sys.argv[1]
out_path = sys.argv[2]

popdict = dict()

with open(fam_path, 'r') as fin:
    with open(out_path, 'w') as fout:
        for line in fin:
            row = line.split()

            # get the existing pop and smaple codes
            long_pop, long_code = row[0:2]

            if long_pop not in popdict:
                # make a new abbreviation
                p = re.split('[^A-Za-z]', long_pop.strip(' _'))

                if len(p) > 1:
                    short_pop = ''.join([w[0] for w in p])
                else:
                    short_pop = long_pop[0:3] if len(long_pop) > 4 else long_pop

                short_pop = short_pop.upper()

                # avoid any duplicate codes
                count = 1
                while short_pop in popdict.values():
                    count += 1
                    short_pop = short_pop + str(count)


                popdict[long_pop] = short_pop

            # get the short names
            row[0] = popdict[long_pop]
            row[1] = long_code.replace(long_pop, short_pop)

            fout.write(' '.join(row) + '\n')

for key in popdict:
    if key != popdict[key]:
        print key, '\t', popdict[key]