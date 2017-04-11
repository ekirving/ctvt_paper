#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, os, os.path, glob, re, shutil, sys

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


class PlinkIndepPairwise(PrioritisedTask):
    """
    Produce a list of SNPs with high discriminating power, by filtering out sites under linkage disequilibrium.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkFilterPops(self.group, self.dataset)

    def output(self):
        extensions = ['in', 'out', 'log']
        return [luigi.LocalTarget("bed/{0}.{1}.prune.{2}".format(self.group, self.dataset, ext)) for ext in extensions]

    def run(self):

        # calculate the prune list (prune.in / prune.out)
        log = run_cmd(["plink",
                       "--dog",
                       "--indep-pairwise", 50, 10, 0.5,  # accept R^2 coefficient of up to 0.5
                       "--bfile", "bed/{0}.{1}".format(self.group, self.dataset),
                       "--out", "bed/{0}.{1}".format(self.group, self.dataset)])

        # write the log file
        with open(self.output()[2].path, 'w') as fout:
            fout.write(log)


class PlinkPruneBed(PrioritisedTask):
    """
    Prune the merged group BED file using the prune.in list from PlinkIndepPairwise()
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkIndepPairwise(self.group, self.dataset)

    def output(self):
        extensions = ['bed', 'bim', 'fam']
        return [luigi.LocalTarget("bed/{0}.{1}.pruned.{2}".format(self.group, self.dataset, ext)) for ext in extensions]

    def run(self):

        # apply the prune list
        run_cmd(["plink",
                 "--dog",
                 "--make-bed",
                 "--extract", self.input()[0].path,
                 "--bfile", "bed/{0}.{1}".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.pruned".format(self.group, self.dataset)])


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


class PlinkBedToFreq(PrioritisedTask):
    """
    Convert a BED file into a minor allele frequency report, needed for input into Treemix.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkFilterGenoByPops(self.group, self.dataset)

    def output(self):
        return luigi.LocalTarget("bed/{0}.{1}.geno.frq.strat.gz".format(self.group, self.dataset))

    def run(self):

        run_cmd(["plink",
                 "--dog",
                 "--freq", "gz",  # make a gzipped MAF report
                 "--family",      # group by population
                 "--bfile", "bed/{0}.{1}.geno".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno".format(self.group, self.dataset)])


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
        return [luigi.LocalTarget("smartpca/{0}.{1}.{2}".format(self.group, self.dataset, ext)) for ext in ['pca.evec', 'eval']]

    def run(self):

        # tell smartpca which pops to use for calculating the eigenvectors, and by inference, which to project
        pops = [pop for pop in GROUPS[self.dataset][self.group] if pop not in self.projectpops]

        # save the pop list
        poplist = 'smartpca/{0}.{1}.poplist'.format(self.group, self.dataset)
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
        ]

        # smartpca needs the params to be defined in a .par file
        parfile = "smartpca/{0}.{1}.par".format(self.group, self.dataset)

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
        for pc1, pc2 in [(1, 2), (3, 4), (5, 6)]:
            yield luigi.LocalTarget("pdf/{0}.{1}.PCA.{2}.{3}.pdf".format(self.group, self.dataset, pc1, pc2))

    def run(self):

        # calculate the percentage of variance explained by each PC, by dividing each eigenvalue by the sum total
        sum = run_cmd(["awk '{ sum+=$1 } END {print sum}' " + self.input()[1].path], shell=True)
        pve = run_cmd(["awk '{print $1/" + sum.strip() +"}' " + self.input()[1].path], shell=True)

        # save the percentage variance
        with open("smartpca/{0}.{1}.pve".format(self.group, self.dataset), 'w') as fout:
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
            with open("smartpca/{0}.{1}.{2}.pca".format(self.group, self.dataset, suffix), 'w') as fout:
                fout.write(data)

        # plot the first 6 components
        for pc1, pc2 in [(1, 2), (3, 4), (5, 6)]:

            # generate both labeled and unlabeled PDFs
            pdfs = {
                0: "pdf/{0}.{1}.PCA.{2}.{3}.pdf".format(self.group, self.dataset, pc1, pc2),
                1: "pdf/{0}.{1}.PCA.{2}.{3}.labeled.pdf".format(self.group, self.dataset, pc1, pc2)
            }

            for labeled, pdf_path in pdfs.iteritems():
                # generate a PDF of the PCA plot
                run_cmd(["Rscript",
                         "rscript/smartpca-plot.R",
                         "smartpca/{0}.{1}.calc.pca".format(self.group, self.dataset),  # pca data, used for calculating the eigenvectors
                         "smartpca/{0}.{1}.proj.pca".format(self.group, self.dataset),  # projected pca data
                         "smartpca/{0}.{1}.pve".format(self.group, self.dataset),       # pve data (% variance)
                         pdf_path,                                                       # location to save the pdf file
                         pc1,                                                            # component num for x-axis
                         pc2,                                                            # component num for y-axis
                         labeled])                                                       # show point labels (0/1)


class TreemixPlinkFreq(PrioritisedTask):
    """
    Convert a Plink MAF frequency file into Treemix format
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkBedToFreq(self.group, self.dataset)

    def output(self):
        return luigi.LocalTarget("treemix/{0}.{1}.geno.frq.gz".format(self.group, self.dataset))

    def run(self):

        # convert the file
        run_cmd(["python",
                 "plink2treemix.py",
                 self.input().path,
                 self.output().path])


class TreemixM(PrioritisedTask):
    """
    Run Treemix with the given number of migration events.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)
    k = luigi.IntParameter(default=1000)

    def requires(self):
        return TreemixPlinkFreq(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("treemix/{0}.{1}.geno.k{2}.m{3}.{4}".format(self.group, self.dataset, self.k, self.m, ext))
                    for ext in ['cov.gz', 'covse.gz', 'edges.gz', 'llik', 'modelcov.gz', 'treeout.gz', 'vertices.gz']]

    def run(self):

        # run treemix
        run_cmd(["treemix",
                 "-i", self.input().path,   # the input file
                 "-root", OUTGROUP_POP[self.dataset],
                 "-k", self.k,              # group together "k" SNPs to account for linkage disequilibrium
                 "-m", self.m,              # build the ML graph with "m" migration events
                 "-o", trim_ext(self.output()[0].path, 2)])


class TreemixPlotM(PrioritisedTask):
    """
    Plot the Treemix output for the given number of migration events.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)
    k = luigi.IntParameter(default=1000)

    def requires(self):
        return TreemixM(self.group, self.dataset, self.m, self.k)

    def output(self):
        return luigi.LocalTarget("pdf/{0}.{1}.treemix.k{2}.m{3}.pdf".format(self.group, self.dataset, self.k, self.m))

    def run(self):

        # compose an ordered population list, with colors for the node labels
        poplist = "treemix/{0}.{1}.poplist".format(self.group, self.dataset)
        with open(poplist, 'w') as fout:
            for pop in GROUPS[self.dataset][self.group]:
                fout.write("{}\t{}\n".format(pop, COLOURS.get(POPULATIONS.get(pop), DEFAULT_COLOUR)))

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
    k = luigi.IntParameter(default=1000)

    def requires(self):
        return TreemixM(self.group, self.dataset, self.m, self.k)

    def output(self):
        return luigi.LocalTarget("qpgraph/{0}.{1}.treemix.k{2}.m{3}.graph".format(self.group, self.dataset, self.k, self.m))

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
    k = luigi.IntParameter(default=1000)

    def requires(self):
        return [TreemixToQPGraph(self.group, self.dataset, self.m, self.k)]

    def output(self):
        return [luigi.LocalTarget("qpgraph/{0}.{1}.treemix.k{2}.m{3}.{4}".format(self.group, self.dataset, self.k, self.m, ext))
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
    k = luigi.IntParameter(default=1000)

    def requires(self):
        return QPGraph(self.group, self.dataset, self.m, self.k)

    def output(self):
        return luigi.LocalTarget("pdf/{0}.{1}.qpgraph.k{2}.m{3}.pdf".format(self.group, self.dataset, self.k, self.m))

    def run(self):

        # conver to postscript format
        ps = run_cmd(["dot", "-Tpdf", self.input()[1].path])

        # save the ps file
        with self.output().open('w') as psfile:
            psfile.write(ps)


class CTVTPipeline(luigi.WrapperTask):
    """
    Run the main CTVC pipeline
    """

    def requires(self):

        for dataset in GROUPS:

            for group in GROUPS[dataset]:

                yield SmartPCAPlot(group, dataset)

                if group not in NO_OUTGROUPS:

                    for m in range(0, TREEMIX_MAX_M + 1):

                        yield TreemixPlotM(group, dataset, m)
                        yield QPGraphPlot(group, dataset, m)


class CTVTCustomPipeline(luigi.WrapperTask):
    """
    Run the specific elements of the CTVC pipeline
    """

    def requires(self):

        dataset = 'merged_map'

        yield SmartPCAPlot('all-pops', dataset)
        yield SmartPCAPlot('all-no-out', dataset)
        yield SmartPCAPlot('dog-ctvt', dataset)
        yield SmartPCAPlot('dog-ctvt', dataset, ['DPC', 'CTVT'])



if __name__ == '__main__':
    luigi.run()