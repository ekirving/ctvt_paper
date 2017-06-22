#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing, copy
from collections import OrderedDict

from luigi import configuration
from psutil import virtual_memory

# populations
ANCIENT_POPS = ['DPC']

ALL_POPS = ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID', 'DQA',
            'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN', 'COY', 'TAI', 'WAM', 'WAS', 'WEU', 'WME', 'OUT']

ALL_DOGS = ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID', 'DQA',
            'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN']

SNP_ARRAY_POPS = ['BAS', 'BAS2', 'DNA', 'BEA', 'BOX', 'DAE', 'DEU', 'DGS', 'DLB', 'AM', 'DAL', 'DGL', 'DHU', 'DMA',
                  'DSL', 'EUR', 'FS', 'GSD', 'SAM', 'SH', 'AED', 'APBT', 'AST', 'CD', 'CLD', 'CBR', 'CHI', 'DME',
                  'DPU', 'NEW', 'NSDTR', 'PIO', 'VDB', 'VDB2', 'VDC', 'VDCR', 'VDDR', 'VDH', 'VDP', 'VDPA', 'VDPC',
                  'VDPI', 'VDPL', 'VDPP', 'VDPR', 'VDUA', 'XOL', 'DID', 'DQA', 'VDIC', 'VDID', 'VDIH', 'VDIM',
                  'VDIO', 'CSP', 'CC', 'DCH', 'DTI', 'DTM', 'DVN', 'NGSD', 'VDIB', 'VDIJ', 'VDPNGEH', 'VDPNGPM',
                  'VDVCB', 'VDVHG', 'VDVLS', 'VDVLC', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME',
                  'TAI', 'OUT']

GRAPH_POPS1 = ['COY', 'WAM', 'WEU', 'DEU', 'DVN', 'DPC', 'DMA']
GRAPH_POPS2 = ['OUT', 'COY', 'WAM', 'WEU', 'DEU', 'DVN', 'DPC', 'DMA']
GRAPH_POPS3 = ['OUT', 'COY', 'WAM', 'WEU', 'TAI', 'DEU', 'DVN', 'DPC', 'DMA']

# groups of populations for running analyses
GROUPS = {

    'merged_map': {

        # all the populations (except the Taimyr)
        'all-pops': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                     'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME', 'OUT'],

        # all the populations, without the outgroup
        'all-no-out': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                       'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME'],

        # dogs + CTVT (no outgroup)
        'dog-ctvt': ['BAS', 'DNA', 'DAE', 'DEU', 'DGS', 'DLB', 'DAL', 'DGL', 'DHU', 'DMA', 'DSL', 'DME', 'DPU', 'DID',
                     'DQA', 'DCH', 'DTI', 'DTM', 'DVN', 'DPC', 'CTVT'],
    },

    'merged_v1':        {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v1.random': {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v1_hq':     {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v1_hq2':    {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v1_TV':     {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v1_TV_hq':  {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},

    'merged_v2':        {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v2.random': {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v2_hq':     {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v2_hq2':    {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v2_TV':     {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v2_TV_hq':  {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},

    # qpGraph analysis
    'merged_v2_laurent':    {'graph-pops1': GRAPH_POPS1,
                             'graph-pops2': GRAPH_POPS2,
                             'graph-pops3': GRAPH_POPS3},

    'merged_v2_TV_laurent': {'graph-pops1': GRAPH_POPS1,
                             'graph-pops2': GRAPH_POPS2,
                             'graph-pops3': GRAPH_POPS3},

    'merged_v3':        {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},
    'merged_v3_TV':     {'all-pops': ALL_POPS, 'dog-ctvt' : ALL_DOGS},


    'merged_SNParray':           {'all-pops': SNP_ARRAY_POPS},
    'merged_SNParray_v1':        {'all-pops': SNP_ARRAY_POPS},
    'merged_SNParray_v1_noCTVT': {'all-pops': SNP_ARRAY_POPS},

}

# add the 'all-pops' group to the new analysis dataset
GROUPS['merged_map_hq']  = {'all-pops': list(GROUPS['merged_map']['all-pops'])}
GROUPS['merged_map_hq2'] = {'all-pops': list(GROUPS['merged_map']['all-pops'])}

# added the Taimyr (Ancient Wolf)
GROUPS['merged_map_Taimyr'] = copy.deepcopy(GROUPS['merged_map'])
for pop in GROUPS['merged_map_Taimyr'] :
    GROUPS['merged_map_Taimyr'][pop].append('Taimyr')

# the population and sample to use for rooting the NJ tree
OUTGROUP_POP = {group: 'OUT' for group in GROUPS}
OUTGROUP_SAMPLE = {group: 'AndeanFox' for group in GROUPS}

# sepcial cases for qpGraph
OUTGROUP_POP['graph-pops1'] = 'COY'
OUTGROUP_SAMPLE['graph-pops1'] = 'C_Cal'

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
    ('AM', 'Arctic Dogs'),
    ('DAL', 'Arctic Dogs'),
    ('DGL', 'Arctic Dogs'),
    ('DHU', 'Arctic Dogs'),
    ('DMA', 'Arctic Dogs'),
    ('DSL', 'Arctic Dogs'),
    ('EUR', 'Arctic Dogs'),
    ('FS', 'Arctic Dogs'),
    ('GSD', 'Arctic Dogs'),
    ('SAM', 'Arctic Dogs'),
    ('SH', 'Arctic Dogs'),
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
    ('DPC', 'Pre-contact Dogs'),
    ('CTVT', 'CTVT'),
    ('DIN', 'Dingo'),
    ('COY', 'Coyote'),
    ('WAM', 'American Wolf'),
    ('WAS', 'Eurasian Wolf'),
    ('WEU', 'Eurasian Wolf'),
    ('WME', 'Eurasian Wolf'),
    ('TAI', 'Ancient Wolf'),
    ('OUT', 'Outgroup')
])

COLOURS = {
    'African Dogs':       '#a6cee3',
    'European Dogs':      '#33a02c',
    'Arctic Dogs':        '#fb9a99',
    'American Dogs':      '#cab2d6',
    'Asian Dogs':         '#fdbf6f',
    'East Asian Dogs':    '#e31a1c',
    'Dingo':              '#ff7f00',
    'Pre-contact Dogs': '#6a3d9a',
    'CTVT':               '#b2df8a',
    'American Wolf':      '#b15928',
    'Eurasian Wolf':      '#003c30',
    'Ancient Wolf':       '#4d4d4d',
    'Coyote':             '#1f78b4',
    'Outgroup':           '#ae017e',
}
DEFAULT_COLOUR = '#e7298a'

# the maximum number of migration events in Treemix
TREEMIX_MAX_M = 5

# number of SNPs to group for LD
TREEMIX_K = 1000

# what level should Treemix group by, pops OR samples?
GROUP_BY_POPS = 'grp-pops'
GROUP_BY_SMPL = 'grp-smpl'

# the maximum number of ancestral populatons to run admixture for
ADMIXTURE_MAX_K = 5

# the number of bootstrap replicates to run
ADMIXTURE_BOOTSTRAP = 0  # TODO put this back to 100

# the species flag for plink telling it how many chromosomes to expect
PLINK_TAXA = '--dog'

# no single worker should use more than 30% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.3)

# dynamically set the resource constraints
config = configuration.get_config()
config.set('resources', 'cpu-cores', str(multiprocessing.cpu_count()))
config.set('resources', 'ram-gb', str(int(virtual_memory().total * 1e-9)))
