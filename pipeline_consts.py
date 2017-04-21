#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing, copy
from collections import OrderedDict

# the population and sample to use for rooting the NJ tree
OUTGROUP_POP = {
    'merged_map': 'OUT',
    'merged_map_hq': 'OUT',
    'merged_map_hq2': 'OUT',
    'merged_map_Taimyr': 'OUT',
    'merged_SNParray': 'OUT'
}

OUTGROUP_SAMPLE = {
    'merged_map': 'AndeanFox',
    'merged_map_hq': 'AndeanFox',
    'merged_map_hq2': 'AndeanFox',
    'merged_map_Taimyr': 'AndeanFox',
    'merged_SNParray': 'AndeanFox'
}

# populations
ANCIENT_POPS = ['DPC']

# groups of populations for running analyses
GROUPS = {

    'merged_map': {

        # all the populations
        'all-pops': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                     'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME', 'OUT'],

        # all the populations, without the outgroup
        'all-no-out': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                       'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME'],

        # dogs + ctvc (no outgroup)
        'dog-ctvt': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                     'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT'],

        # TODO remove when done testing
        'test-pops' : ['OUT', 'CTVT', 'DPC', 'WAM']
    },

    'merged_SNParray': {

        'all-pops' : ['BAS', 'BAS2', 'DNA', 'BEA', 'BOX', 'DAE', 'DEU', 'DGS', 'DLB', 'AM', 'DAL', 'DGL', 'DHU', 'DMA',
                      'DSL', 'EUR', 'FS', 'GSD', 'SAM', 'SH', 'AED', 'APBT', 'AST', 'CD', 'CLD', 'CBR', 'CHI', 'DME',
                      'DPU', 'NEW', 'NSDTR', 'PIO', 'VDB', 'VDB2', 'VDC', 'VDCR', 'VDDR', 'VDH', 'VDP', 'VDPA', 'VDPC',
                      'VDPI', 'VDPL', 'VDPP', 'VDPR', 'VDUA', 'XOL', 'DID', 'DQA', 'VDIC', 'VDID', 'VDIH', 'VDIM',
                      'VDIO', 'CSP', 'CC', 'DCH', 'DTI', 'DTM', 'DVN', 'NGSD', 'VDIB', 'VDIJ', 'VDPNGEH', 'VDPNGPM',
                      'VDVCB', 'VDVHG', 'VDVLS', 'VDVLC', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME',
                      'TAI', 'OUT']
    },

}

# add the 'all-pops' group to the new analysis dataset
GROUPS['merged_map_hq']  = {'all-pops': list(GROUPS['merged_map']['all-pops'])}
GROUPS['merged_map_hq2'] = {'all-pops': list(GROUPS['merged_map']['all-pops'])}

# added the Taimyr (Ancient Wolf)
GROUPS['merged_map_Taimyr'] = copy.deepcopy(GROUPS['merged_map'])
for pop in GROUPS['merged_map_Taimyr'] :
    GROUPS['merged_map_Taimyr'][pop].append('Taimyr')


NO_OUTGROUPS = ['all-no-out', 'dog-ctvt']

POPULATIONS = OrderedDict([
    ('BAS', 'African Dogs'),
    ('BAS2', 'African Dogs'),
    ('DNA', 'African Dogs'),
    ('BEA', 'European Dogs'),
    ('BOX', 'European Dogs'),
    ('DAE', 'European Dogs'),
    ('DEU', 'European Dogs'),
    ('DGS', 'European Dogs'),
    ('DLB', 'European Dogs'),
    ('AM', 'Northern Dogs'),
    ('DAL', 'Northern Dogs'),
    ('DGL', 'Northern Dogs'),
    ('DHU', 'Northern Dogs'),
    ('DMA', 'Northern Dogs'),
    ('DSL', 'Northern Dogs'),
    ('EUR', 'Northern Dogs'),
    ('FS', 'Northern Dogs'),
    ('GSD', 'Northern Dogs'),
    ('SAM', 'Northern Dogs'),
    ('SH', 'Northern Dogs'),
    ('AED', 'American Dogs'),
    ('APBT', 'American Dogs'),
    ('AST', 'American Dogs'),
    ('CD', 'American Dogs'),
    ('CLD', 'American Dogs'),
    ('CBR', 'American Dogs'),
    ('CHI', 'American Dogs'),
    ('DME', 'American Dogs'),
    ('DPU', 'American Dogs'),
    ('NEW', 'American Dogs'),
    ('NSDTR', 'American Dogs'),
    ('PIO', 'American Dogs'),
    ('VDB', 'American Dogs'),
    ('VDB2', 'American Dogs'),
    ('VDC', 'American Dogs'),
    ('VDCR', 'American Dogs'),
    ('VDDR', 'American Dogs'),
    ('VDH', 'American Dogs'),
    ('VDP', 'American Dogs'),
    ('VDPA', 'American Dogs'),
    ('VDPC', 'American Dogs'),
    ('VDPI', 'American Dogs'),
    ('VDPL', 'American Dogs'),
    ('VDPP', 'American Dogs'),
    ('VDPR', 'American Dogs'),
    ('VDUA', 'American Dogs'),
    ('XOL', 'American Dogs'),
    ('DID', 'Asian Dogs'),
    ('DQA', 'Asian Dogs'),
    ('VDIC', 'Asian Dogs'),
    ('VDID', 'Asian Dogs'),
    ('VDIH', 'Asian Dogs'),
    ('VDIM', 'Asian Dogs'),
    ('VDIO', 'Asian Dogs'),
    ('CSP', 'East Asian Dogs'),
    ('CC', 'East Asian Dogs'),
    ('DCH', 'East Asian Dogs'),
    ('DTI', 'East Asian Dogs'),
    ('DTM', 'East Asian Dogs'),
    ('DVN', 'East Asian Dogs'),
    ('NGSD', 'East Asian Dogs'),
    ('VDIB', 'East Asian Dogs'),
    ('VDIJ', 'East Asian Dogs'),
    ('VDPNGEH', 'East Asian Dogs'),
    ('VDPNGPM', 'East Asian Dogs'),
    ('VDVCB', 'East Asian Dogs'),
    ('VDVHG', 'East Asian Dogs'),
    ('VDVLS', 'East Asian Dogs'),
    ('VDVLC', 'East Asian Dogs'),
    ('DPC', 'Pre-Colombian Dogs'),
    ('CTVT', 'CTVT'),
    ('DIN', 'Dingo'),
    ('COY', 'Coyotes'),
    ('WAM', 'American Wolf'),
    ('WAS', 'Eurasian Wolf'),
    ('WEU', 'Eurasian Wolf'),
    ('WME', 'Eurasian Wolf'),
    ('TAI', 'Ancient Wolf'),
    ('OUT', 'Outgroup')
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
    'Outgroup':           '#4d4d4d',
    'American Wolf':      '#b15928',
    'Ancient Wolf':       '#4d4d4d',
    'Eurasian Wolf':      '#003c30',
}
DEFAULT_COLOUR = '#e7298a'

# the maximum number of migration events in Treemix
TREEMIX_MAX_M = 10

# number of SNPs to group for LD
TREEMIX_K = 1000

# what level should Treemix group by, pops OR samples?
GROUP_BY_POPS = 'grp-pops'
GROUP_BY_SMPL = 'grp-smpl'

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)
