#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

# import my custom modules
from pipeline_utils import *

if len(sys.argv) != 3:
    print "Error: Missing PED files"
    quit()

ped_path = sys.argv[1]
out_path = sys.argv[2]

# make a random allele call for each site (i.e. pretend everything is homo)
parse_ped(ped_path + '.ped',
          out_path + '.ped')

# copy the map file
copyfile(ped_path + '.map',
         out_path + '.map')