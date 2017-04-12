#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
from collections import OrderedDict

# the population and sample to use for rooting the NJ tree
OUTGROUP_POP = {
    'merged_map': 'OUT',
    'merged_map_hq': 'OUT'
}

OUTGROUP_SAMPLE = {
    'merged_map': 'AndeanFox',
    'merged_map_hq': 'AndeanFox'
}

# populations
ANCIENT_POPS = ['DPC']

# groups of populations for running analyses
GROUPS = {

    'merged_map': {

        # all the populations
        'all-pops': ['BAS', 'COY', 'CTVT', 'DAE', 'DAL', 'DCH', 'DEU', 'DGL', 'DGS', 'DHU', 'DID', 'DIN', 'DLB', 'DMA',
                     'DME', 'DNA', 'DPC', 'DPU', 'DQA', 'DSL', 'DTI', 'DTM', 'DVN', 'OUT', 'WAM', 'WAS', 'WEU', 'WME'],

        # all the populations, without the outgroup
        'all-no-out': ['BAS', 'COY', 'CTVT', 'DAE', 'DAL', 'DCH', 'DEU', 'DGL', 'DGS', 'DHU', 'DID', 'DIN', 'DLB', 'DMA',
                       'DME', 'DNA', 'DPC', 'DPU', 'DQA', 'DSL', 'DTI', 'DTM', 'DVN', 'WAM', 'WAS', 'WEU', 'WME'],

        # dogs + ctvc (no outgroup)
        'dog-ctvt': ['BAS', 'CTVT', 'DAE', 'DAL', 'DCH', 'DEU', 'DGL', 'DGS', 'DHU', 'DID', 'DLB', 'DMA', 'DME', 'DNA',
                     'DPC', 'DPU', 'DQA', 'DSL', 'DTI', 'DTM', 'DVN'],

        'test-pops' : ['OUT', 'CTVT', 'DPC', 'WAM']
    },
}

# add the 'all-pops' group to the new analysis dataset
GROUPS['merged_map_hq']  = {'all-pops': GROUPS['merged_map']['all-pops']}
GROUPS['merged_map_hq2'] = {'all-pops': GROUPS['merged_map']['all-pops']}

NO_OUTGROUPS = ['all-no-out', 'dog-ctvt']

POPULATIONS = OrderedDict([
    ('BAS', 'African Dogs'),
    ('COY', 'Coyotes'),
    ('CTVT', 'CTVT'),
    ('DAE', 'European Dogs'),
    ('DAL', 'Northern Dogs'),
    ('DCH', 'East Asian Dogs'),
    ('DEU', 'European Dogs'),
    ('DGL', 'Northern Dogs'),
    ('DGS', 'European Dogs'),
    ('DHU', 'Northern Dogs'),
    ('DID', 'Asian Dogs'),
    ('DIN', 'Dingo'),
    ('DLB', 'European Dogs'),
    ('DMA', 'Northern Dogs'),
    ('DME', 'American Dogs'),
    ('DNA', 'African Dogs'),
    ('DPC', 'Pre-Colombian Dogs'),
    ('DPU', 'American Dogs'),
    ('DQA', 'Asian Dogs'),
    ('DSL', 'Northern Dogs'),
    ('DTI', 'East Asian Dogs'),
    ('DTM', 'East Asian Dogs'),
    ('DVN', 'East Asian Dogs'),
    ('OUT', 'Outgroup'),
    ('WAM', 'American Wolf'),
    ('WAS', 'Eurasian Wolf'),
    ('WEU', 'Eurasian Wolf'),
    ('WME', 'Eurasian Wolf'),
])

COLOURS = {
    'African Dogs':       '#a6cee3',
    'Coyotes':            '#1f78b4',
    'CTVT':               '#b2df8a',
    'European Dogs':      '#33a02c',
    'Northern Dogs':      '#fb9a99',
    'East Asian Dogs':    '#e31a1c',
    'Asian Dogs':         '#fdbf6f',
    'Dingo':              '#ff7f00',
    'American Dogs':      '#cab2d6',
    'Pre-Colombian Dogs': '#6a3d9a',
    'Outgroup':           '#ffff99',
    'American Wolf':      '#b15928',
    'Eurasian Wolf':      '#003c30',
}
DEFAULT_COLOUR = '#e7298a'

# the samtools flag for BAM file compression
DEFAULT_COMPRESSION = 6

# the minimum phred scaled genotype quality (30 = 99.9%)
MIN_BASE_QUAL = 30
MIN_MAPPING_QUAL = 30

# number of bases to soft clip
SOFT_CLIP_DIST = 5

# column delimiter for PED files
PED_COL_DELIM = ' '

# status codes for PED files
PED_UNKNOWN = '0'
PED_MISSING_PHENO = '-9'
PED_MISSING_GENO = '0'

BIM_COL_CHROM = 0
BIM_COL_POS = 3

# the maximum number of ancestral populatons to run admixture for
ADMIXTURE_MAX_K = 5

# the number of bootstrap replicates to run
ADMIXTURE_BOOTSTRAP = 0  # TODO put this back to 100

# the maximum number of migration events in Treemix
TREEMIX_MAX_M = 10

# number of SNPs to group for LD
TREEMIX_K = 1000

# what level should Treemix group by, pops OR samples?
GROUP_BY_POPS = 'grp-pops'
GROUP_BY_SMPL = 'grp-smpl'

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)

# location of software tools
SAMTOOLS = "/usr/local/bin/samtools1.3"
TRIMMOMATIC = "/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar"
TRIMMOMATIC_ADAPTERS = "/usr/local/Trimmomatic-0.36/adapters/TruSeq3-SE.fa"
