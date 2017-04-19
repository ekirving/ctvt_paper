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

        'all-pops' : ['BAS', 'Basenji', 'DNA', 'Beagle', 'Boxer', 'DAE', 'DEU', 'DGS', 'DLB', 'Alaskan_Malamute', 'DAL',
                      'DGL', 'DHU', 'DMA', 'DSL', 'Eurasier', 'Finnish_Spitz', 'Greenland_Sledge_Dog', 'Samoyed',
                      'Siberian_Husky', 'American_Eskimo_Dog', 'American_Pit_Bull_Terrier', 'American_Staffordshire_Terrier',
                      'Carolina_Dog', 'Catahoula_Leopard_Dog', 'Chesapeake_Bay_Retriever', 'Chihuahua', 'DME', 'DPU',
                      'Newfoundland', 'Nova_Scotia_Duck_Tolling_Retriever', 'Peruvian_Inca_Orchid', 'Village_Dog_Belize',
                      'Village_Dog_Brazil', 'Village_Dog_Colombia', 'Village_Dog_Costa_Rica', 'Village_Dog_Dominican_Republic',
                      'Village_Dog_Honduras', 'Village_Dog_Panama', 'Village_Dog_Peru-Arequipa', 'Village_Dog_Peru-Cusco',
                      'Village_Dog_Peru-Ica', 'Village_Dog_Peru-Loreto', 'Village_Dog_Peru-Puno', 'Village_Dog_Puerto_Rico',
                      'Village_Dog_US-Alaska', 'Xoloitzcuintli', 'DID', 'DQA', 'Village_Dog_India-Chennai',
                      'Village_Dog_India-Dehli', 'Village_Dog_India-Hazaribagh', 'Village_Dog_India-Mumbai',
                      'Village_Dog_India-Orissa', 'Chinese_Shar-pei', 'Chow_Chow', 'DCH', 'DTI', 'DTM', 'DVN',
                      'New_Guinea_Singing_Dog', 'Village_Dog_Indonesia-Borneo', 'Village_Dog_Indonesia-Jakarta',
                      'Village_Dog_Papua_New_Guinea-East_Highlands_', 'Village_Dog_Papua_New_Guinea-Port_Moresby',
                      'Village_Dog_Vietnam-Cao_Bang', 'Village_Dog_Vietnam-Ha_Giang', 'Village_Dog_Vietnam-Lang_Son',
                      'Village_Dog_Vietnam-Lao_Cai', 'DPC', 'CTVT', 'DIN', 'COY', 'WAM', 'WAS', 'WEU', 'WME', 'Taimyr', 'OUT']
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
    'Outgroup':           '#4d4d4d',
    'American Wolf':      '#b15928',
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