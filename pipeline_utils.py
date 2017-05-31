#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, subprocess, datetime, hashlib, os, random, itertools

# import all the constants
from pipeline_consts import *

def parse_ped(ped_infile, ped_outfile):
    """
    Randomly select allele at heterozygous site in ped file and convert to homozygous
    """
    with open(ped_infile, 'r') as fin:
        with open(ped_outfile, 'w') as fout:

            for line in fin:
                row = line.split()
                info = row[0:6]
                genotype_rand = row[6:]
                it = iter(genotype_rand)
                genotype_it = zip(it,it)
                for n,i in enumerate(genotype_it):
                    rand_allele=random.choice(i)
                    genotype_it[n]=rand_allele,rand_allele
                genotype_t = list(itertools.chain(*genotype_it))
                new_line = info+genotype_t

                fout.write(" ".join(new_line)+"\n")


def run_cmd(cmd, returnout=True, shell=False, pwd='./', verbose=True):
    """
    Executes the given command in a system subprocess

    :param cmd: The system command to run (list|string)
    :param shell: Use the native shell
    :return: The stdout stream
    """
    # subprocess only accepts strings
    cmd = [str(args) for args in cmd]

    # print the command
    if verbose:
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


def get_metapops(group, dataset, metalist):
    """
    Get all the populations belonging to the given meta groups.
    """
    return [pop for pop in GROUPS[dataset][group] if POPULATIONS[pop] in metalist]


def pprint_qpgraph(dot_file, pdf_file):
    """
    Pretty print a qpGraph dot file using the graphviz library.
    """
    import graphviz
    import re

    # extract the body contents from the dot file
    with open(dot_file, 'rU') as fin:
        body = re.search("{(.+)}", fin.read(), re.DOTALL).group(1).strip().split("\n")

    # make a new direct graph
    dot = graphviz.Digraph(filename=trim_ext(pdf_file), body=body, format='pdf')

    # remove the messy graph label
    dot.attr('graph', label='')

    # set Node defaults
    dot.node_attr['shape'] = 'point'
    dot.node_attr['fontname'] = 'arial'
    dot.node_attr['fontsize'] = '11'

    # set Edge defaults
    dot.edge_attr['arrowhead'] = 'vee'
    dot.edge_attr['fontcolor'] = '#838b8b'  # grey
    dot.edge_attr['fontname'] = 'arial'
    dot.edge_attr['fontsize'] = '11'

    nodes = []

    # extract the leaf nodes from the body of the graph
    for line in dot:
        match = re.search("^ +([a-z]+) +\[", line, re.IGNORECASE)
        if match:
            nodes.append(match.group(1))

    # set leaf node attributes
    for node in nodes:
        colour = COLOURS.get(POPULATIONS.get(node), DEFAULT_COLOUR)
        dot.node(node, shape='ellipse', color=colour, fontcolor=colour)

    # render the graph (basename.pdf)
    dot.render(cleanup=True)


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
