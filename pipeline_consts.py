#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
from collections import OrderedDict

# the population and sample to use for rooting the NJ tree
OUTGROUP_POP = {}

OUTGROUP_SAMPLE = {}

# populations
ANCIENT_POPS = ['DPC']

# groups of populations for running analyses
GROUPS = {

    'ascertain1': {

        # all the populations
        'all-pops': {'BAS', 'COY', 'CTVT', 'DAE', 'DAL', 'DCH', 'DEU', 'DGL', 'DGS', 'DHU', 'DID', 'DIN', 'DLB', 'DMA',
                    'DME', 'DNA', 'DPC', 'DPU', 'DQA', 'DSL', 'DTI', 'DTM', 'DVN', 'OUT', 'WAM', 'WAS', 'WEU', 'WME'},
    }
}

NO_OUTGROUPS = []

POPULATIONS = OrderedDict([

])

COLOURS = {

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

# no single worker should use more than 50% of the available cores
MAX_CPU_CORES = int(multiprocessing.cpu_count() * 0.5)

# location of software tools
SAMTOOLS = "/usr/local/bin/samtools1.3"
TRIMMOMATIC = "/usr/local/Trimmomatic-0.36/trimmomatic-0.36.jar"
TRIMMOMATIC_ADAPTERS = "/usr/local/Trimmomatic-0.36/adapters/TruSeq3-SE.fa"
