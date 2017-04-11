#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, subprocess, datetime, hashlib, os, pysam, random

from shutil import copyfile

# import all the constants
from pipeline_consts import *

ped = open(sys.argv[1],"r")
#ped = open("toy.ped","r")
outfile = open("test.ped","w")
hom = [('A','A'),('C','C'),('G','G'),('T','T')]

def parse_ped(ped_infile, ped_outfile):
    """
    randomly select allele at heterozygous site in ped file and convert to homozygous
    """
    hom = [('A', 'A'), ('C', 'C'), ('G', 'G'), ('T', 'T')]

    for line in ped_infile:
        row = line.split()
        info = row[0:6]
        genotype_rand = row[6:]
        it = iter(genotype_rand)
        genotype_it = zip(it,it)
        for n,i in enumerate(genotype_it):
                if i not in hom:
                    rand_allele=random.choice(i)
                    genotype_it[n]=rand_allele,rand_allele
        genotype_t = list(itertools.chain(*genotype_it))
        new_line = info+genotype_t
        ped_outfile.write(" ".join(new_line)+"\n")
    ped_outfile.close()


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
