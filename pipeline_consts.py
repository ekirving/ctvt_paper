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
    ('Basenji', 'African Dogs'),
    ('DNA', 'African Dogs'),
    ('Beagle', 'European Dogs'),
    ('Boxer', 'European Dogs'),
    ('DAE', 'European Dogs'),
    ('DEU', 'European Dogs'),
    ('DGS', 'European Dogs'),
    ('DLB', 'European Dogs'),
    ('Alaskan_Malamute', 'Northern Dogs'),
    ('DAL', 'Northern Dogs'),
    ('DGL', 'Northern Dogs'),
    ('DHU', 'Northern Dogs'),
    ('DMA', 'Northern Dogs'),
    ('DSL', 'Northern Dogs'),
    ('Eurasier', 'Northern Dogs'),
    ('Finnish_Spitz', 'Northern Dogs'),
    ('Greenland_Sledge_Dog', 'Northern Dogs'),
    ('Samoyed', 'Northern Dogs'),
    ('Siberian_Husky', 'Northern Dogs'),
    ('American_Eskimo_Dog', 'American Dogs'),
    ('American_Pit_Bull_Terrier', 'American Dogs'),
    ('American_Staffordshire_Terrier', 'American Dogs'),
    ('Carolina_Dog', 'American Dogs'),
    ('Catahoula_Leopard_Dog', 'American Dogs'),
    ('Chesapeake_Bay_Retriever', 'American Dogs'),
    ('Chihuahua', 'American Dogs'),
    ('DME', 'American Dogs'),
    ('DPU', 'American Dogs'),
    ('Newfoundland', 'American Dogs'),
    ('Nova_Scotia_Duck_Tolling_Retriever', 'American Dogs'),
    ('Peruvian_Inca_Orchid', 'American Dogs'),
    ('Village_Dog_Belize', 'American Dogs'),
    ('Village_Dog_Brazil', 'American Dogs'),
    ('Village_Dog_Colombia', 'American Dogs'),
    ('Village_Dog_Costa_Rica', 'American Dogs'),
    ('Village_Dog_Dominican_Republic', 'American Dogs'),
    ('Village_Dog_Honduras', 'American Dogs'),
    ('Village_Dog_Panama', 'American Dogs'),
    ('Village_Dog_Peru-Arequipa', 'American Dogs'),
    ('Village_Dog_Peru-Cusco', 'American Dogs'),
    ('Village_Dog_Peru-Ica', 'American Dogs'),
    ('Village_Dog_Peru-Loreto', 'American Dogs'),
    ('Village_Dog_Peru-Puno', 'American Dogs'),
    ('Village_Dog_Puerto_Rico', 'American Dogs'),
    ('Village_Dog_US-Alaska', 'American Dogs'),
    ('Xoloitzcuintli', 'American Dogs'),
    ('DID', 'Asian Dogs'),
    ('DQA', 'Asian Dogs'),
    ('Village_Dog_India-Chennai', 'Asian Dogs'),
    ('Village_Dog_India-Dehli', 'Asian Dogs'),
    ('Village_Dog_India-Hazaribagh', 'Asian Dogs'),
    ('Village_Dog_India-Mumbai', 'Asian Dogs'),
    ('Village_Dog_India-Orissa', 'Asian Dogs'),
    ('Chinese_Shar-pei', 'East Asian Dogs'),
    ('Chow_Chow', 'East Asian Dogs'),
    ('DCH', 'East Asian Dogs'),
    ('DTI', 'East Asian Dogs'),
    ('DTM', 'East Asian Dogs'),
    ('DVN', 'East Asian Dogs'),
    ('New_Guinea_Singing_Dog', 'East Asian Dogs'),
    ('Village_Dog_Indonesia-Borneo', 'East Asian Dogs'),
    ('Village_Dog_Indonesia-Jakarta', 'East Asian Dogs'),
    ('Village_Dog_Papua_New_Guinea-East_Highlands_', 'East Asian Dogs'),
    ('Village_Dog_Papua_New_Guinea-Port_Moresby', 'East Asian Dogs'),
    ('Village_Dog_Vietnam-Cao_Bang', 'East Asian Dogs'),
    ('Village_Dog_Vietnam-Ha_Giang', 'East Asian Dogs'),
    ('Village_Dog_Vietnam-Lang_Son', 'East Asian Dogs'),
    ('Village_Dog_Vietnam-Lao_Cai', 'East Asian Dogs'),
    ('DPC', 'Pre-Colombian Dogs'),
    ('CTVT', 'CTVT'),
    ('DIN', 'Dingo'),
    ('COY', 'Coyotes'),
    ('WAM', 'American Wolf'),
    ('WAS', 'Eurasian Wolf'),
    ('WEU', 'Eurasian Wolf'),
    ('WME', 'Eurasian Wolf'),
    ('Taimyr', 'Ancient Wolf'),
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