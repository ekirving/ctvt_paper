#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, os, os.path, glob, re, shutil, sys, random, itertools

# import my custom modules
from pipeline_utils import *

from collections import defaultdict


class AscertainedBed(luigi.ExternalTask):
    """
    External task that checks that starting bed file has been supplied.
    """
    dataset = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}".format(self.dataset, ext)) for ext in ['bed', 'bim', 'fam']]


class PlinkFilterPops(PrioritisedTask):
    """
    Filter out unwanted populations from a BED file.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return AscertainedBed(self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.{2}".format(self.group, self.dataset, ext)) for ext in ['bed', 'bim', 'fam']]

    def run(self):
        poplist = "bed/{0}.{1}.poplist".format(self.group, self.dataset)

        # write the population IDs to disk so we can filter out any unwanted "families"
        with open(poplist, 'w') as par:
            par.write("\n".join(GROUPS[self.dataset][self.group]))

        # apply the prune list
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 "--keep-fam", poplist,
                 "--bfile", trim_ext(self.input()[0].path),
                 "--out", "bed/{0}.{1}".format(self.group, self.dataset)])


class PlinkExtractPop(PrioritisedTask):
    """
    Extract a single population from a larger BED file
    """
    dataset = luigi.Parameter()
    population = luigi.Parameter()

    def requires(self):
        return AscertainedBed(self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.{2}".format(self.dataset, self.population, ext)) for ext in ['bed', 'bim', 'fam']]

    def run(self):
        poplist = "bed/{0}.{1}.poplist".format(self.dataset, self.population)

        # write the population IDs to disk so we can filter out any unwanted "families"
        with open(poplist, 'w') as par:
            par.write(self.population)

        # apply the prune list
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 "--geno", "0.999",  # drop sites with no coverage, or rather, where less than 0.1% is missing
                 "--keep-fam", poplist,
                 "--bfile", trim_ext(self.input()[0].path),
                 "--out", "bed/{0}.{1}".format(self.dataset, self.population)])


class PlinkFilterGenoByPops(PrioritisedTask):
    """
    Filter a merged BED file, such that all sites are called at least once in all the sub-populations of interest.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        yield PlinkFilterPops(self.group, self.dataset)

        for population in ANCIENT_POPS:
            yield PlinkExtractPop(self.dataset, population)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/{0}.{1}.geno.{2}".format(self.group, self.dataset, ext)) for ext in
                extensions]

    def run(self):

        vars = []

        # get the list of all the variant IDs which have coverage in each ancient population
        for pop in ANCIENT_POPS:
            bimfile = "bed/{0}.{1}.bim".format(self.dataset, pop)
            vars.append([line.split()[1] for line in open(bimfile).readlines()])

        # get the intersection of those variant lists (i.e. only sites covered in all populations)
        unqvars = set(vars[0]).intersection(*vars)

        # save the variant list
        varlist = 'bed/{0}.{1}.varlist'.format(self.group, self.dataset)
        with open(varlist, 'w') as fout:
            fout.write("\n".join(unqvars))

        # filter the big BED file using that list of variants
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 "--extract", varlist,
                 "--bfile", "bed/{0}.{1}".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno".format(self.group, self.dataset)])


class RandomPedAllele(PrioritisedTask):
    """
    Convert bed to ped. Convert heterozygous sites in ped file to homozygous choosing alleles at random.
    Convert ped back to bed
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkFilterGenoByPops(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.geno.random.{2}".format(self.group, self.dataset, ext)) for ext in
                    ['bed', 'bim', 'fam']]

    def run(self):

        # conver to PED so it's easier to parse the data
        run_cmd(["plink",
                 "--dog",
                 "--recode",
                 "--bfile", "bed/{0}.{1}.geno".format(self.group, self.dataset),
                 "--out", "ped/{0}.{1}.geno".format(self.group, self.dataset),])

        # make a random allele call for each site (i.e. pretend everything is homo)
        parse_ped("ped/{0}.{1}.geno.ped".format(self.group, self.dataset),
                  "ped/{0}.{1}.geno.random.ped".format(self.group, self.dataset))

        # copy the map file
        copyfile("ped/{0}.{1}.geno.map".format(self.group, self.dataset),
                 "ped/{0}.{1}.geno.random.map".format(self.group, self.dataset))

        # convert the random called file back to BED
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 "--file", "ped/{0}.{1}.geno.random".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno.random".format(self.group, self.dataset)])


class PlinkHighGeno(PrioritisedTask):
    """
    Filter all sites without a gentyping rate of 95%
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return RandomPedAllele(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.geno.random.hq.{2}".format(self.group, self.dataset, ext)) for ext in
                    ['bed', 'bim', 'fam']]

    def run(self):

        # conver to PED so it's easier to parse the data
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 # "--maf", "0.05",
                 "--geno", "0.05",
                 "--bfile", "bed/{0}.{1}.geno.random".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno.random.hq".format(self.group, self.dataset)])


class PlinkBedToFreq(PrioritisedTask):
    """
    Convert a BED file into a minor allele frequency report, needed for input into Treemix.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkHighGeno(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.geno.random.{2}.{3}".format(self.group, self.dataset, self.groupby, ext))
                    for ext in ['frq.strat.gz', 'log']]

    def run(self):

        # get the paths for the input files
        bed_file, bim_file, fam_file = [file.path for file in self.input()]

        # to group by samples, we need to reassign them to their own families
        if self.groupby == GROUP_BY_SMPL:

            # replace family with sample code
            fam = run_cmd(["awk '{$1=$2}$0' " + fam_file], shell=True)

            # make a new fam file
            fam_file = insert_suffix(fam_file, GROUP_BY_SMPL)
            with open(fam_file, 'w') as fout:
                fout.write(fam)

        run_cmd(["plink",
                 "--dog",
                 "--freq", "gz",     # make a gzipped MAF report
                 "--family",         # group by population
                 "--bed", bed_file,
                 "--bim", bim_file,
                 "--fam", fam_file,
                 "--out", "bed/{0}.{1}.geno.random.{2}".format(self.group, self.dataset, self.groupby)])

class SmartPCA(PrioritisedTask):
    """
    Calcualte the eigenvectors for the given population
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    projectpops = luigi.ListParameter(default=ANCIENT_POPS)

    def requires(self):
        return PlinkFilterPops(self.group, self.dataset)

    def output(self):
        prj = "-".join(self.projectpops)
        return [luigi.LocalTarget("smartpca/{0}.{1}.prj-{2}.{3}".format(self.group, self.dataset, prj, ext)) for ext in ['pca.evec', 'eval']]

    def run(self):

        # tell smartpca which pops to use for calculating the eigenvectors, and by inference, which to project
        pops = [pop for pop in GROUPS[self.dataset][self.group] if pop not in self.projectpops]

        prj = "-".join(self.projectpops)

        # save the pop list
        poplist = 'smartpca/{0}.{1}.prj-{2}.poplist'.format(self.group, self.dataset, prj)
        with open(poplist, 'w') as fout:
            fout.write("\n".join(pops))

        # N.B. smartpca requires an invalid verison of the .fam file, where the phenotype column has been replaced with
        # the population name. If you don't do this then none of the populations in poplistname will be found. Also, if
        # you happen to be using qtmode:YES (because otherwise a value of -9 throws errors!) then it will fail to find
        # the populations, even if you've make the invalid version of the .fam file that poplistname wants. Grrrr... >:(
        fam = run_cmd(["awk '$6=$1' " + self.input()[2].path], shell=True)

        # save the fam file
        famfile = insert_suffix(self.input()[2].path, "smartpca")
        with open(famfile, 'w') as fout:
            fout.write(fam)

        # compose the config settings for smartpca
        config = [
            "genotypename:   {0}".format(self.input()[0].path),   # .bed
            "snpname:        {0}".format(self.input()[1].path),   # .bim
            "indivname:      {0}".format(famfile),                # invalid .fam
            "evecoutname:    {0}".format(self.output()[0].path),  # .pca.evec
            "evaloutname:    {0}".format(self.output()[1].path),  # .eval
            "numthreads:     {0}".format(MAX_CPU_CORES),          # number of threads to use
            "poplistname:    {0}".format(poplist),                # the list of pops to calculate the eigenvectors from
            "lsqproject:     YES",                                # use least squares projection, best for missing data
            "numoutlieriter: 0",                                  # don't exclude outliers
            "inbreed:        YES",                                # compute inbreeding stats
        ]

        # smartpca needs the params to be defined in a .par file
        parfile = "smartpca/{0}.{1}.prj-{2}.par".format(self.group, self.dataset, prj)

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # calcualt the PCA and project the ancient samples ontop
        run_cmd(["smartpca",
                 "-p", parfile])


class SmartPCAPlot(PrioritisedTask):
    """
    Use ggplot to plot the PCA
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    projectpops = luigi.ListParameter(default=ANCIENT_POPS)

    def requires(self):
        return SmartPCA(self.group, self.dataset, self.projectpops)

    def output(self):
        prj = "-".join(self.projectpops)
        for pc1, pc2 in [(1, 2), (3, 4), (5, 6)]:
            yield luigi.LocalTarget("pdf/{0}.{1}.prj-{2}.PCA.{3}.{4}.pdf".format(self.group, self.dataset, prj, pc1, pc2))

    def run(self):
        prj = "-".join(self.projectpops)

        # calculate the percentage of variance explained by each PC, by dividing each eigenvalue by the sum total
        sum = run_cmd(["awk '{ sum+=$1 } END {print sum}' " + self.input()[1].path], shell=True)
        pve = run_cmd(["awk '{print $1/" + sum.strip() +"}' " + self.input()[1].path], shell=True)

        # save the percentage variance
        with open("smartpca/{0}.{1}.prj-{2}.pve".format(self.group, self.dataset, prj), 'w') as fout:
            fout.write(pve)

        # copy the population name into the first column, skip the header row (to make things easier for the Rscript)
        awk = "awk 'NR>1 {print $NF $0}' " + self.input()[0].path

        # make a regex to match only the ancient pops
        projectpops = "|".join(["^{} ".format(pop) for pop in self.projectpops])

        # split the data into calculated and projected
        calc = run_cmd([awk + ' | grep -Ev "{}"'.format(projectpops)], shell=True)
        proj = run_cmd([awk + ' | grep -E  "{}"'.format(projectpops)], shell=True)

        # save the pca data
        for suffix, data in [('calc', calc), ('proj',proj)]:
            with open("smartpca/{0}.{1}.prj-{2}.{3}.pca".format(self.group, self.dataset, prj, suffix), 'w') as fout:
                fout.write(data)

        # plot the first 6 components
        for pc1, pc2 in [(1, 2), (3, 4), (5, 6)]:

            # generate both labeled and unlabeled PDFs
            pdfs = {
                0: "pdf/{0}.{1}.prj-{2}.PCA.{3}.{4}.pdf".format(self.group, self.dataset, prj, pc1, pc2),
                1: "pdf/{0}.{1}.prj-{2}.PCA.{3}.{4}.labeled.pdf".format(self.group, self.dataset, prj, pc1, pc2)
            }

            for labeled, pdf_path in pdfs.iteritems():
                # generate a PDF of the PCA plot
                run_cmd(["Rscript",
                         "rscript/smartpca-plot.R",
                         "smartpca/{0}.{1}.prj-{2}.calc.pca".format(self.group, self.dataset, prj),  # pca data, used for calculating the eigenvectors
                         "smartpca/{0}.{1}.prj-{2}.proj.pca".format(self.group, self.dataset, prj),  # projected pca data
                         "smartpca/{0}.{1}.prj-{2}.pve".format(self.group, self.dataset, prj),       # pve data (% variance)
                         pdf_path,                                                       # location to save the pdf file
                         pc1,                                                            # component num for x-axis
                         pc2,                                                            # component num for y-axis
                         labeled])                                                       # show point labels (0/1)


class NeighborJoiningTree(PrioritisedTask):
    """
    Create a neighbor joining phylogenetic tree from a pruned BED file
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    treetype = luigi.Parameter(default='phylogram')

    def requires(self):
        return PlinkFilterGenoByPops(self.group, self.dataset)

    def output(self):
            return [luigi.LocalTarget("njtree/{0}.{1}.geno.data".format(self.group, self.dataset)),
                    luigi.LocalTarget("njtree/{0}.{1}.geno.tree".format(self.group, self.dataset)),
                    luigi.LocalTarget("pdf/{0}.{1}.njtree.pdf".format(self.group, self.dataset))]

    def run(self):

        # TODO what about bootstrapping?
        # make the distance matrix
        run_cmd(["plink",
                 "--dog",
                 "--distance", "square", "1-ibs",
                 "--bfile", "bed/{0}.{1}.geno".format(self.group, self.dataset),
                 "--out", "njtree/{0}.{1}.geno".format(self.group, self.dataset)])

        fam_file = "bed/{0}.{1}.geno.fam".format(self.group, self.dataset)

        # use awk to extract the sample names
        awk = "awk '{print $2}' " + fam_file

        # transpose them into a row
        head = run_cmd([awk + " | xargs"], returnout=True, shell=True)

        # use awk to extract the sample and population names
        awk = "awk '{print $1\"\\t\"$2}' " + fam_file

        # add the samples names as a column to the mdist data
        data = run_cmd([awk + " | paste - njtree/{0}.{1}.geno.mdist".format(self.group, self.dataset)], shell=True)

        # save the labeled file
        with self.output()[0].open('w') as fout:
            fout.write("Code\tSample\t" + head)
            fout.write(data)

        # generate a tree from the labeled data
        run_cmd(["Rscript",
                 "rscript/plot-phylo-tree.R",
                 self.output()[0].path,
                 self.treetype,
                 OUTGROUP_SAMPLE[self.dataset],
                 self.output()[1].path,
                 self.output()[2].path])


class TreemixPlinkFreq(PrioritisedTask):
    """
    Convert a Plink MAF frequency file into Treemix format
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkBedToFreq(self.group, self.dataset, self.groupby)

    def output(self):
        return luigi.LocalTarget("treemix/{0}.{1}.geno.random.{2}.frq.gz".format(self.group, self.dataset, self.groupby))

    def run(self):

        # convert the file
        run_cmd(["python",
                 "plink2treemix.py",
                 self.input()[0].path,
                 self.output().path])


class TreemixM(PrioritisedTask):
    """
    Run Treemix with the given number of migration events.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    groupby = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return TreemixPlinkFreq(self.group, self.dataset, self.groupby)

    def output(self):
        return [luigi.LocalTarget("treemix/{0}.{1}.geno.random.{2}.m{3}.{4}".format(self.group, self.dataset, self.groupby, self.m, ext))
                    for ext in ['cov.gz', 'covse.gz', 'edges.gz', 'llik', 'modelcov.gz', 'treeout.gz', 'vertices.gz']]

    def run(self):

        if self.groupby == GROUP_BY_POPS:
            outgroup = OUTGROUP_POP[self.dataset]
        else:
            outgroup = OUTGROUP_SAMPLE[self.dataset]

        # run treemix
        run_cmd(["treemix",
                 "-i", self.input().path,   # the input file
                 "-root", outgroup,
                 "-k", TREEMIX_K,           # group together "k" SNPs to account for linkage disequilibrium
                 "-m", self.m,              # build the ML graph with "m" migration events
                 "-o", trim_ext(self.output()[0].path, 2)])


class TreemixPlotM(PrioritisedTask):
    """
    Plot the Treemix output for the given number of migration events.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    groupby = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return TreemixM(self.group, self.dataset, self.groupby, self.m)

    def output(self):
        return luigi.LocalTarget("pdf/{0}.{1}.treemix.geno.random.{2}.m{3}.pdf".format(self.group, self.dataset, self.groupby, self.m))

    def run(self):

        # compose an ordered population list, with colors for the node labels
        poplist = "treemix/{0}.{1}.geno.random.{2}.poplist".format(self.group, self.dataset, self.groupby)

        with open(poplist, 'w') as fout:
            if self.groupby == GROUP_BY_POPS:
                # output the populations
                for pop in GROUPS[self.dataset][self.group]:
                    colour = COLOURS.get(POPULATIONS.get(pop), DEFAULT_COLOUR)
                    fout.write("{}\t{}\n".format(pop, colour))
            else:
                # fetch the sample names from the fam file
                fam_file = "bed/{0}.{1}.geno.random.fam".format(self.group, self.dataset)
                fam = run_cmd(["awk '{print $2\" \"$1}' " + fam_file], shell=True)
                samples = dict(line.split() for line in fam.splitlines())

                # output the samples
                for sample in samples:
                    colour = COLOURS.get(POPULATIONS.get(samples[sample]), DEFAULT_COLOUR)
                    fout.write("{}\t{}\n".format(sample, colour))

        # plot the treemix tree
        run_cmd(["Rscript",
                 "rscript/treemix-plot.R",
                 trim_ext(self.input()[0].path, 2),
                 poplist,
                 self.output().path])


class TreemixToQPGraph(PrioritisedTask):
    """
    Convert a Treemix model into qpGraph format
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return TreemixM(self.group, self.dataset, GROUP_BY_POPS, self.m)

    def output(self):
        return luigi.LocalTarget("qpgraph/{0}.{1}.treemix.m{2}.graph".format(self.group, self.dataset, self.m))

    def run(self):

        # the path to treemix treeout file
        treeout = self.input()[5].path

        # perform the conversion
        graph = run_cmd(["python", "tree2qpg.py", treeout])

        with self.output().open('w') as fout:
            fout.write(graph)


class QPGraph(PrioritisedTask):
    """
    Run qpGraph
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return [TreemixToQPGraph(self.group, self.dataset, self.m)]

    def output(self):
        return [luigi.LocalTarget("qpgraph/{0}.{1}.treemix.m{2}.{3}".format(self.group, self.dataset, self.m, ext))
                    for ext in ['par', 'dot', 'log']]

    def run(self):

        # NB admixtools requires a non-standard fam file format
        fam = run_cmd(["awk '$6=$1' bed/{0}.{1}.geno.fam".format(self.group, self.dataset)], shell=True)

        # save the fam file
        famfile = "bed/{0}.{1}.qpgraph.fam".format(self.group, self.dataset)
        with open(famfile, 'w') as fout:
            fout.write(fam)

        # TODO comment these options
        # compose the config settings for qpGraph
        config = [
            "genotypename:  bed/{0}.{1}.geno.bed".format(self.group, self.dataset),
            "snpname:       bed/{0}.{1}.geno.bim".format(self.group, self.dataset),
            "indivname:     bed/{0}.{1}.qpgraph.fam".format(self.group, self.dataset),
            "outpop:        {}".format(OUTGROUP_POP[self.dataset]),
            "blgsize:       0.05",
            "lsqmode:       YES",
            "diag:          .0001",
            "hires:         YES",
            "initmix:       1000",
            "precision:     .0001",
            "zthresh:       3.0",
            "terse:         NO",
            "useallsnps:    NO",
        ]

        # qpGraph needs the params to be defined in a .par file
        parfile = self.output()[0].path

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # run qpGraph
        log = run_cmd(["qpGraph",
                       "-p", parfile,
                       "-g", self.input()[0].path,    # graph file
                       "-d", self.output()[1].path])  # dot file

        # save the log file
        with self.output()[2].open('w') as logfile:
            logfile.write(log)


class QPGraphPlot(PrioritisedTask):
    """
    Plot the output from qpGraph
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return QPGraph(self.group, self.dataset, self.m)

    def output(self):
        return luigi.LocalTarget("pdf/{0}.{1}.qpgraph.m{2}.pdf".format(self.group, self.dataset, self.m))

    def run(self):

        # conver to postscript format
        ps = run_cmd(["dot", "-Tpdf", self.input()[1].path])

        # save the ps file
        with self.output().open('w') as psfile:
            psfile.write(ps)


class ConvertfBedToEigenstrat(PrioritisedTask):
    """
    Convert a BED file into Eigenstrat format, for use by admixtools
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkFilterPops(self.group, self.dataset)

    def output(self):
        extensions = ['par', 'eigenstratgeno', 'snp', 'ind', 'log']
        return [luigi.LocalTarget("eigenstrat/{0}.{1}.{2}".format(self.group, self.dataset, ext)) for ext in extensions]

    def run(self):

        # NB admixtools requires a non-standard fam file format
        # also... we want to pretent that all our samples are populations!
        fam = run_cmd(["awk '$6=$2' " + self.input()[2].path], shell=True)

        # save the fam file
        famfile = insert_suffix(self.input()[2].path, "convertf")
        with open(famfile, 'w') as fout:
            fout.write(fam)

        # compose the config settings for convertf
        config = [
            "genotypename:    {}".format(self.input()[0].path),
            "snpname:         {}".format(self.input()[1].path),
            "indivname:       {}".format(famfile),
            "outputformat:    EIGENSTRAT",
            "genotypeoutname: {}".format(self.output()[1].path),
            "snpoutname:      {}".format(self.output()[2].path),
            "indivoutname:    {}".format(self.output()[3].path),
            "familynames:     NO",
            "pordercheck:     NO"
        ]

        # convertf needs the params to be defined in a .par file
        parfile = self.output()[0].path

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # run convertf
        log = run_cmd(["convertf",
                       "-p", parfile])

        # save the log file
        with self.output()[4].open('w') as logfile:
            logfile.write(log)


class QPDstat(PrioritisedTask):
    """
    Run qpDstat to test all four-population admixture models.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    blgsize = luigi.Parameter()

    # only run one at a time
    resources = {'qpdstat': 1}

    def requires(self):
        return ConvertfBedToEigenstrat(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("qpdstat/{0}.{1}.blgsize-{2}.{3}".format(self.group, self.dataset, self.blgsize, ext))
                    for ext in ['par', 'log', 'poplist']]

    def run(self):

        # get the out group pop
        outpop = OUTGROUP_POP[self.dataset]

        famfile = "bed/{0}.{1}.convertf.fam".format(self.group, self.dataset)

        # split the samples based on outgroup
        out_samples = run_cmd(["grep    '" + outpop + "' " + famfile + " | awk '{print $2}'"], shell=True).splitlines()
        all_samples = run_cmd(["grep -v '" + outpop + "' " + famfile + " | awk '{print $2}'"], shell=True).splitlines()

        perms = itertools.permutations(all_samples, 3)

        # write the list of 4-way tests
        with self.output()[2].open('w') as fout:

            for out in out_samples:
                for pop1, pop2, pop3 in perms:
                    fout.write(" ".join([out, pop1, pop2, pop3]) + "\n")

        # compose the config settings
        config = [
            "genotypename: {}".format(self.input()[1].path),
            "snpname:      {}".format(self.input()[2].path),
            "indivname:    {}".format(self.input()[3].path),
            "popfilename:  {}".format(self.output()[2].path),  # Program will run the method for all listed quadrapules
            "blgsize:      {}".format(self.blgsize)
            # "f4mode:   YES",                                   # f4 statistics not D-stats are computed
        ]

        # the params to be defined in a .par file
        parfile = self.output()[0].path

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # run qpDstat
        log = run_cmd(["qpDstat", "-p", parfile])

        # save the log file
        with self.output()[1].open('w') as logfile:
            logfile.write(log)


class CTVTPipeline(luigi.WrapperTask):
    """
    Run the main CTVT pipeline
    """

    def requires(self):

        for dataset in GROUPS:

            for group in GROUPS[dataset]:

                yield SmartPCAPlot(group, dataset)

                for blgsize in [0.5, 1, 2]:
                    yield QPDstat(group, dataset, blgsize)

                if group not in NO_OUTGROUPS:

                    for m in range(0, TREEMIX_MAX_M + 1):

                        yield TreemixPlotM(group, dataset, GROUP_BY_POPS, m)
                        yield TreemixPlotM(group, dataset, GROUP_BY_SMPL, m)
                        yield QPGraphPlot(group, dataset, m)


class CTVTCustomPipeline(luigi.WrapperTask):
    """
    Run the specific elements of the CTVT pipeline
    """

    def requires(self):

        # all the data

        for dataset in ['merged_map', 'merged_map_Taimyr', 'merged_SNParray']:

            yield SmartPCAPlot('all-pops', dataset)
            yield NeighborJoiningTree('all-pops', dataset)

        for dataset in ['merged_map', 'merged_map_Taimyr']:

            yield SmartPCAPlot('all-no-out', dataset)
            yield SmartPCAPlot('dog-ctvt', dataset)
            yield SmartPCAPlot('dog-ctvt', dataset, ['DPC', 'CTVT'])

            # don't attempt this with merged_SNParray as there are >258e6 combinations
            for blgsize in [0.5, 1, 2]:
                yield QPDstat('all-pops', dataset, blgsize)

        # only the high quality ancient samples
        for dataset in ['merged_map_hq', 'merged_map_hq2', 'merged_SNParray']:
            for m in range(0, 5):
                yield TreemixPlotM('all-pops', dataset, GROUP_BY_POPS, m)
                if dataset != 'merged_SNParray':
                    yield TreemixPlotM('all-pops', dataset, GROUP_BY_SMPL, m)


class CTVTCustomPipelineV2(luigi.WrapperTask):
    """
    Run the specific elements of the CTVT pipeline
    """

    def requires(self):

        for dataset in ['merged_v1', 'merged_v2']:

            yield SmartPCAPlot('all-pops', dataset, ANCIENT_POPS)
            yield SmartPCAPlot('all-pops', dataset + '.random', ANCIENT_POPS)

            yield NeighborJoiningTree('all-pops', dataset)



if __name__ == '__main__':
    luigi.run()
