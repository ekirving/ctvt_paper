#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, subprocess, datetime, hashlib, os, pysam, random

from shutil import copyfile

# import all the constants
from pipeline_consts import *


def run_cmd(cmd, returnout=True, shell=False, pwd='./'):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # TODO remove when done testing
    print cmd

    # has the command so we can match the logs together
    m = hashlib.md5()
    m.update(str(cmd))

    # log the command
    with open(pwd + 'log/luigi.cmd.log', 'a+') as fout:
        fout.write("{time} {hash}# {actn}\n".format(time=str(datetime.datetime.now()),
                                                    hash=m.hexdigest(),
                                                    actn=" ".join(cmd)))

    # run the command
    proc = subprocess.Popen(cmd,
                            shell=shell,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # fetch the output and error
    (stdout, stderr) = proc.communicate()

    # bail if something went wrong
    if proc.returncode:
        raise Exception(stderr)

    if returnout:
        return stdout
    else:
        # TODO send output to a log file
        # with open('', 'w') as fout:
        pass


def unzip_file(gzip):
    """
    Unzip a gzipped file, using multi-threading when available

    :param gzip: The path to the gzip file
    :return: The uncompressed data stream
    """
    try:
        # use unpigz for multithreaded unzipping (if installed)
        return run_cmd(["unpigz",
                        "-c",                 # output to stdout
                        "-p", MAX_CPU_CORES,  # use the maximum cores
                        gzip])                # gzip file

    except OSError as e:
        # if it's not installed
        if e.errno == os.errno.ENOENT:

            # use single-threaded gunzip
            return run_cmd(["gunzip",
                            "-c",             # output to stdout
                            gzip])            # gzip file
        else:
            # escalate the exception
            raise e


def curl_download(url, filename):
    """
    Downloads a remote url to a local file path using cURL
    """

    # download the file
    run_cmd(["curl",
             "-s",                  # download silently
             "--output", filename,  # output path
             url])                  # from this url


def generate_ped_file(population, sites):
    """
    Generates a biallelic PED file by aggregating base calls, in a given list of sites, found in multiple BAM files.

    Each BAM file is converted to a pileup, quality filtered, then randomly selected (if multiple bases remain).

    :param population: The population code
    :param sites: Dictionary of (chrom, pos) to include in the PED file
    :return:
    """

    # process all the BAM files
    for sample in ANCIENT_POPS[population]:

        # open the BAM file for reading
        with pysam.AlignmentFile('bam/{}.sorted.bam'.format(sample), 'rb') as bamfile:

            # iterate over the whole file, in pileup mode (this is faster than making >700k random access requests)
            for pileupcolumn in bamfile.pileup():

                # get the chrom and pos of the current site in the BAM file
                chrom, pos = pileupcolumn.reference_name, pileupcolumn.reference_pos

                # skip any sites not found in the SNP array
                if (chrom, pos) not in sites:
                    # TODO error logging
                    continue

                quality_bases = []
                clipped_bases = []

                # iterate over all the reads for this site
                for pileupread in pileupcolumn.pileups:

                    # skip alignments that don't have a base at this site (e.g. InDels)
                    if not pileupread.is_del and not pileupread.is_refskip:

                        # skip bad alignments
                        if pileupread.alignment.mapping_quality < MIN_MAPPING_QUAL:
                            # TODO error logging
                            continue

                        # skip bad base calls
                        if pileupread.alignment.query_qualities[pileupread.query_position] < MIN_BASE_QUAL:
                            # TODO error logging
                            continue

                        # get the aligned base for this read
                        base = pileupread.alignment.query_sequence[pileupread.query_position]

                        # get the overall length of the read
                        read_length = len(pileupread.alignment.query_sequence)

                        # soft clip bases near the edges of the read
                        if pileupread.query_position > SOFT_CLIP_DIST or \
                            pileupread.query_position < (read_length - SOFT_CLIP_DIST):

                            quality_bases.append(base)
                        else:
                            clipped_bases.append(base)

                # randomly select a base, giving preference to quality bases over soft clipped
                if len(quality_bases) > 0:
                    sites[(chrom, pos)][sample] = random.choice(quality_bases)

                elif len(clipped_bases) > 0:
                    # TODO error logging
                    sites[(chrom, pos)][sample] = random.choice(clipped_bases)

                else:
                    # TODO error logging
                    pass

    # PED can only handle biallelic sites
    for site in sites.keys():

        # so set polyallelic sites to unknown
        if len(set(sites[site].values())) > 2:
            # TODO error logging
            sites[site] = {}

    ped = ''

    # process one sample at a time
    for sample in ANCIENT_POPS[population]:
        row = list()

        # Family ID
        row.append(population)

        # Individual ID
        row.append(sample)

        # Paternal ID
        row.append(PED_UNKNOWN)

        # Maternal ID
        row.append(PED_UNKNOWN)

        # Sex
        row.append(PED_UNKNOWN)

        # Phenotype
        row.append(PED_MISSING_PHENO)

        # Genotypes
        for site in sites:

            if sample in sites[site]:
                # add the genotype twice, because our coverage is too low to call het sites
                for _ in range(2):
                    row.append(sites[site][sample])
            else:
                # TODO error logging
                for _ in range(2):
                    row.append(PED_MISSING_GENO)

        row.append('\n')

        # print the row to the file
        ped += PED_COL_DELIM.join(row)

    return ped


def trim_ext(fullpath, n=1):
    return ('.').join(fullpath.split('.')[:-n])


def trim_path_ext(fullpath):
    return trim_ext(os.path.basename(fullpath))

def insert_suffix(fullpath, suffix):
    splitpath = fullpath.split('.')
    splitpath.insert(-1, suffix)
    return ('.').join(splitpath)


class PrioritisedTask(luigi.Task):
    """
    PrioritisedTask that implements a dynamic priority method
    """
    @property
    def priority(self):

        p = 100

        try:
            # deprioritise large values of k
            p -= self.k
        except AttributeError:
            pass

        try:
            # deprioritise large values of m
            p -= self.m
        except AttributeError:
            pass

        return p
