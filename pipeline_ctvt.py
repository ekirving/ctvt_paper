#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi, os, os.path, glob, re, shutil, sys, random, itertools

# import my custom modules
from pipeline_utils import *
from permute_qpgraph import *

from shutil import copyfile
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
                 PLINK_TAXA,
                 "--make-bed",
                 "--keep-fam", poplist,
                 "--bfile", trim_ext(self.input()[0].path),
                 "--out", "bed/{0}.{1}".format(self.group, self.dataset)])


class PlinkIndepPairwise(PrioritisedTask):
    """
    Produce a list of SNPs with high discriminating power, by filtering out minor allele frequency and sites under
    linkage disequilibrium
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
                       PLINK_TAXA,
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
                 PLINK_TAXA,
                 "--make-bed",
                 "--extract", self.input()[0].path,
                 "--bfile", "bed/{0}.{1}".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.pruned".format(self.group, self.dataset)])


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
                 PLINK_TAXA,
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
                 PLINK_TAXA,
                 "--make-bed",
                 "--extract", varlist,
                 "--bfile", "bed/{0}.{1}".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno".format(self.group, self.dataset)])


class PlinkHighGeno(PrioritisedTask):
    """
    Filter all sites without a genotyping rate of 95%
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return PlinkFilterGenoByPops(self.group, self.dataset)

    def output(self):
        return [luigi.LocalTarget("bed/{0}.{1}.geno.hq.{2}".format(self.group, self.dataset, ext)) for ext in
                    ['bed', 'bim', 'fam']]

    def run(self):

        # conver to PED so it's easier to parse the data
        run_cmd(["plink",
                 PLINK_TAXA,
                 "--make-bed",
                 # "--maf", "0.05",
                 "--geno", "0.05",
                 "--bfile", "bed/{0}.{1}.geno".format(self.group, self.dataset),
                 "--out", "bed/{0}.{1}.geno.hq".format(self.group, self.dataset)])


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
        return [luigi.LocalTarget("bed/{0}.{1}.geno.{2}.{3}".format(self.group, self.dataset, self.groupby, ext))
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
                 PLINK_TAXA,
                 "--freq", "gz",  # make a gzipped MAF report
                 "--family",  # group by population
                 "--bed", bed_file,
                 "--bim", bim_file,
                 "--fam", fam_file,
                 "--out", "bed/{0}.{1}.geno.{2}".format(self.group, self.dataset, self.groupby)])

class SmartPCA(PrioritisedTask):
    """
    Calcualte the eigenvectors for the given population
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    projectpops = luigi.ListParameter(default=ANCIENT_POPS)

    resources = {'cpu-cores': CPU_CORES_MED}

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

        # TODO replace with dependency on convertf
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
            "numthreads:     {0}".format(CPU_CORES_MED),          # number of threads to use
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
    components = luigi.ListParameter(default=PCA_COMPONENTS)

    def requires(self):
        return SmartPCA(self.group, self.dataset, self.projectpops)

    def output(self):
        prj = "-".join(self.projectpops)
        for pc1, pc2 in self.components:
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




class AdmixtureK(PrioritisedTask):
    """
    Run admixture, with K ancestral populations, on the pruned reference data BED file.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    k = luigi.IntParameter()

    resources = {'cpu-cores': CPU_CORES_MED}

    def requires(self):
        return PlinkPruneBed(self.group, self.dataset)

    def output(self):
        extensions = ['P', 'Q', 'log']
        return [luigi.LocalTarget("admix/{0}.{1}.pruned.{2}.{3}".format(self.group, self.dataset, self.k, ext))
                    for ext in extensions]

    def run(self):

        # admixture only outputs to the current directory
        os.chdir('./admix')

        bed_file = "../{0}".format(self.input()[0].path)

        log = run_cmd(["admixture", 
                       "-j{0}".format(CPU_CORES_MED),       # use multi-threading
                       "-B{}".format(ADMIXTURE_BOOTSTRAP),  # the number of bootstrap replicates to run
                       "--cv=10",                           # generate cross-validation estimates
                       bed_file,                            # using the pruned data file
                       self.k],                             # for K ancestral populations
                      pwd='../')

        # restore previous working directory
        os.chdir('..')

        # save the log file
        with self.output()[2].open('w') as fout:
            fout.write(log)


class AdmixtureSortK(PrioritisedTask):
    """
    Sort the admixture Q matrix for the given value of K, such that the column orders maintain optimal relative
    relationships.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    k = luigi.IntParameter()

    def requires(self):
        yield AdmixtureK(self.group, self.dataset, self.k)

        # sorting of k requires that k-1 also be sorted (because we need something to compare against)
        if self.k > 1:
            yield AdmixtureSortK(self.group, self.dataset, self.k - 1)

    def output(self):
        return luigi.LocalTarget("admix/{0}.{1}.pruned.sorted.{2}.Q".format(self.group, self.dataset, self.k))

    def run(self):

        # no need to sort a file with only one k
        if self.k == 1:
            copyfile(self.input()[0][1].path, self.output().path)
            return

        # sort all the matricies
        run_cmd(["Rscript",
                 "rscript/admix-sort-k.R",
                 self.input()[0][1].path,
                 self.input()[1].path,
                 self.output().path])


class AdmixturePlotK(PrioritisedTask):
    """
    Use ggplot to plot the admixture Q stats
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    maxk = luigi.IntParameter()
    k = luigi.IntParameter()

    def requires(self):
        return AdmixtureSortK(self.group, self.dataset, self.k)

    def output(self):
        return [luigi.LocalTarget("admix/{0}.{1}.pruned.sorted.{2}.data".format(self.group, self.dataset, self.k)),
                luigi.LocalTarget("pdf/{0}.{1}.admix.K.{2}.pdf".format(self.group, self.dataset, self.k))]

    def run(self):

        fam = "bed/{0}.{1}.pruned.fam".format(self.group, self.dataset)
        q = "admix/{0}.{1}.pruned.sorted.{2}.Q".format(self.group, self.dataset, self.k)

        # use awk and paste to add population and sample names, needed for the plot
        awk = "awk '{ print $1 \" \" $2 }' " + fam + " | paste - " + q

        # get the admix data
        data = run_cmd([awk], returnout=True, shell=True)

        # parse the data into a dict, so it can be output in a specific order
        datadict = defaultdict(list)
        for line in data.strip().split("\n"):
            datadict[line.split()[0]].append(line)

        # compose the header row
        header = ["Pop{}".format(i) for i in range(1, int(self.k) + 1)]
        header.insert(0, "Sample")
        header.insert(0, "Population")

        # save the labeled file
        with self.output()[0].open('w') as fout:
            # output the header row
            fout.write("\t".join(header)+"\n")

            # output the populations, in the chosen order
            for pop in GROUPS[self.dataset][self.group]:
                for line in datadict[pop]:
                    fout.write(line+"\n")

        # generate a PDF of the admixture stacked column chart
        run_cmd(["Rscript",
                 "rscript/admix-plot-k.R",
                 self.output()[0].path,
                 self.output()[1].path,
                 self.maxk])


class AdmixtureCV(PrioritisedTask):
    """
    Run admixture for the given population, determine the optimal K value, and plot the graphs
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        # set maximum K to be the number of actual populations + 10
        maxk = len(GROUPS[self.dataset][self.group]) + ADMIXTURE_MAX_K

        # run admixture or each population and each value of K
        for k in range(1, maxk + 1):
            yield AdmixturePlotK(self.group, self.dataset, maxk, k)

    def output(self):
        return [luigi.LocalTarget("admix/{0}.{1}.pruned.CV.data".format(self.group, self.dataset)),
                luigi.LocalTarget("pdf/{0}.{1}.admix.CV.pdf".format(self.group, self.dataset))]

    def run(self):

        # use grep to extract the cross-validation scores from all the log files
        cvs = run_cmd(["grep -h CV admix/{0}.{1}.*.log".format(self.group, self.dataset)], returnout=True, shell=True)

        # extract the K value and CV score
        # e.g. "CV error (K=1): 1.20340"
        data = [tuple([re.sub(r'[^\d.]+', "", val) for val in cv.split(":")]) for cv in cvs.splitlines()]

        # get the three lowest CV scores
        # bestfit = sorted(data, key=lambda x: x[1])[0:3]

        # write the scores to a data file
        with self.output()[0].open('w') as fout:
            fout.write("\t".join(["K", "CV"]) + "\n")
            for row in data:
                fout.write("\t".join(str(datum) for datum in row) + "\n")

        # plot the CV values as a line graph
        run_cmd(["Rscript",
                 "rscript/plot-line-graph.R",
                 self.output()[0].path,
                 self.output()[1].path,
                 "Ancestral populations (K)",
                 "Cross-validation Error"])


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
                 PLINK_TAXA,
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

        outgroup = OUTGROUP_SAMPLE[self.group] if self.group in OUTGROUP_SAMPLE else OUTGROUP_SAMPLE[self.dataset]

        # generate a tree from the labeled data
        run_cmd(["Rscript",
                 "rscript/plot-phylo-tree.R",
                 self.output()[0].path,
                 self.treetype,
                 outgroup,
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
        return luigi.LocalTarget("treemix/{0}.{1}.geno.{2}.frq.gz".format(self.group, self.dataset, self.groupby))

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
        return [luigi.LocalTarget("treemix/{0}.{1}.geno.{2}.m{3}.{4}".format(self.group, self.dataset, self.groupby, self.m, ext))
                    for ext in ['cov.gz', 'covse.gz', 'edges.gz', 'llik', 'modelcov.gz', 'treeout.gz', 'vertices.gz']]

    def run(self):

        if self.groupby == GROUP_BY_POPS:
            outgroup = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]
        else:
            outgroup = OUTGROUP_SAMPLE[self.group] if self.group in OUTGROUP_SAMPLE else OUTGROUP_SAMPLE[self.dataset]

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
        return luigi.LocalTarget("pdf/{0}.{1}.treemix.geno.{2}.m{3}.pdf".format(self.group, self.dataset, self.groupby, self.m))

    def run(self):

        # compose an ordered population list, with colors for the node labels
        poplist = "treemix/{0}.{1}.geno.{2}.poplist".format(self.group, self.dataset, self.groupby)

        with open(poplist, 'w') as fout:
            if self.groupby == GROUP_BY_POPS:
                # output the populations
                for pop in GROUPS[self.dataset][self.group]:
                    colour = COLOURS.get(POPULATIONS.get(pop), DEFAULT_COLOUR)
                    fout.write("{}\t{}\n".format(pop, colour))
            else:
                # fetch the sample names from the fam file
                fam_file = "bed/{0}.{1}.geno.fam".format(self.group, self.dataset)
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


class QPGraphTreemix(PrioritisedTask):
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

        outgroup = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]

        # compose the config settings for qpGraph
        config = [
            "genotypename:  bed/{0}.{1}.geno.bed".format(self.group, self.dataset),
            "snpname:       bed/{0}.{1}.geno.bim".format(self.group, self.dataset),
            "indivname:     bed/{0}.{1}.qpgraph.fam".format(self.group, self.dataset),
            "outpop:        {}".format(outgroup),
            # TODO review these defaults
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


class QPGraphTreemixPlot(PrioritisedTask):
    """
    Plot the output from qpGraph
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    m = luigi.IntParameter(default=0)

    def requires(self):
        return QPGraphTreemix(self.group, self.dataset, self.m)

    def output(self):
        return luigi.LocalTarget("pdf/{0}.{1}.qpgraph.m{2}.pdf".format(self.group, self.dataset, self.m))

    def run(self):

        dot_file = self.input()[1].path
        pdf_file = self.output().path

        # pretty print the qpgraph dot file
        pprint_qpgraph(dot_file, pdf_file)


class QPGraphPermute(PrioritisedTask):
    """
    Permute all possible graphs to find a qpGraph model which fits the data
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    exhaustive = luigi.BoolParameter(default=False)

    resources = {'cpu-cores': CPU_CORES_MAX}

    def requires(self):
        return ConvertfBedToEigenstrat(self.group, self.dataset, GROUP_BY_POPS)

    def output(self):
        return [luigi.LocalTarget("qpgraph/{0}.{1}.permute.{2}".format(self.group, self.dataset, ext))
                for ext in ['par', 'log', 'fitted']]

    def run(self):

        # get the populations and the outgroup
        populations = GROUPS[self.dataset][self.group]
        outgroup = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]

        # compose the config settings for qpGraph
        config = [
            "genotypename:  {}".format(self.input()[1].path),
            "snpname:       {}".format(self.input()[2].path),
            "indivname:     {}".format(self.input()[3].path),
            "outpop:        NULL",  # NB. if the outgroup is inbred then all f2 stats are fucked up
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
        par_file = self.output()[0].path

        with open(par_file, 'w') as par:
            par.write("\n".join(config))

        dot_path = 'qpgraph/dot/{0}.permute'.format(self.dataset)
        pdf_path = 'pdf/{0}.{1}.qpg-permute'.format(self.group, self.dataset)

        log_file = self.output()[1].path

        # run qpGraph
        solutions = permute_qpgraph(par_file, log_file, dot_path, pdf_path, populations, outgroup, self.exhaustive,
                                    nthreads=CPU_CORES_MAX)

        # record the hash codes for all the fitting graphs
        fitted_file = self.output()[2].path
        with open(fitted_file, 'w') as fout:
            fout.write("\n".join(solutions))


class QPGraphCluster(PrioritisedTask):
    """
    Perform hierarchical clustering on a set of fitted admixture graphs.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    exhaustive = luigi.BoolParameter(default=False)

    resources = {'cpu-cores': CPU_CORES_HIGH}

    def requires(self):
        return QPGraphPermute(self.group, self.dataset, self.exhaustive)

    def output(self):
        log_file = 'qpgraph/{0}.{1}.cluster.log'.format(self.group, self.dataset)
        csv_file = 'qpgraph/{0}.{1}.cluster.csv'.format(self.group, self.dataset)
        mtx_file = 'qpgraph/{0}.{1}.cluster.npy'.format(self.group, self.dataset)
        pdf_file = 'pdf/{0}.{1}.qpg-cluster.pdf'.format(self.group, self.dataset)

        return [luigi.LocalTarget(path) for path in [log_file, csv_file, mtx_file, pdf_file]]

    def run(self):

        # the prefix to apply to the dot files
        dot_path = 'qpgraph/dot/{0}.permute'.format(self.dataset)

        # get the output file
        log_file, csv_file, mtx_file, pdf_file = [file.path for file in self.output()]

        fitted_file = self.input()[2].path
        with open(fitted_file, 'r') as fin:
            graph_names = fin.read().splitlines()

        cluster_qpgraph(graph_names, dot_path, log_file, pdf_file, csv_file, mtx_file, nthreads=CPU_CORES_HIGH)


class ConvertfBedToEigenstrat(PrioritisedTask):
    """
    Convert a BED file into Eigenstrat format, for use by admixtools
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    groupby = luigi.Parameter()

    def requires(self):
        return PlinkFilterPops(self.group, self.dataset)

    def output(self):
        extensions = ['par', 'eigenstratgeno', 'snp', 'ind', 'log']
        return [luigi.LocalTarget("eigenstrat/{0}.{1}.{2}.{3}".format(self.group, self.dataset, self.groupby, ext)) for ext in extensions]

    def run(self):

        # NB admixtools requires a non-standard fam file format
        if self.groupby == GROUP_BY_POPS:
            fam = run_cmd(["awk '$6=$1' " + self.input()[2].path], shell=True)

        elif self.groupby == GROUP_BY_SMPL:
            # also... we want to pretend that all our samples are populations!
            fam = run_cmd(["awk '$6=$2' " + self.input()[2].path], shell=True)

        # save the fam file
        famfile = insert_suffix(self.input()[2].path, self.groupby)
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

    resources = {'ram-gb': 64}

    def requires(self):
        return ConvertfBedToEigenstrat(self.group, self.dataset, GROUP_BY_SMPL)

    def output(self):
        return [luigi.LocalTarget("qpdstat/{0}.{1}.blgsize-{2}.{3}".format(self.group, self.dataset, self.blgsize, ext))
                    for ext in ['par', 'log', 'poplist']]

    def run(self):

        # get the out group pop
        outpop = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]

        famfile = "bed/{0}.{1}.{2}.fam".format(self.group, self.dataset, GROUP_BY_SMPL)

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


class QP3Pop(PrioritisedTask):
    """
    Run qp3Pop to test all three-population admixture models.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()

    def requires(self):
        return ConvertfBedToEigenstrat(self.group, self.dataset, GROUP_BY_POPS)

    def output(self):
        return [luigi.LocalTarget("qp3pop/{0}.{1}.{2}".format(self.group, self.dataset, ext))
                    for ext in ['par', 'log', 'poplist']]

    def run(self):

        # get the out group pop
        outpop = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]

        # write the list of 3-way tests
        with self.output()[2].open('w') as fout:

            for target in ANCIENT_POPS:
                for testpop in GROUPS[self.dataset][self.group]:
                    if testpop not in [target, outpop]:
                        fout.write(" ".join([outpop, testpop, target]) + "\n")

        # compose the config settings
        config = [
            "genotypename: {}".format(self.input()[1].path),
            "snpname:      {}".format(self.input()[2].path),
            "indivname:    {}".format(self.input()[3].path),
            "popfilename:  {}".format(self.output()[2].path),
            "inbreed:      YES"  # Use if target pop is inbred OR abd crucially if target is psudo-diploid
        ]

        # the params need to be defined in a .par file
        parfile = self.output()[0].path

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # run qp3pop
        log = run_cmd(["qp3Pop", "-p", parfile])

        # save the log file
        with self.output()[1].open('w') as logfile:
            logfile.write(log)


class QPF4ratio(PrioritisedTask):
    """
    Run qpF4ratio to test all F4 ratios.
    """
    group = luigi.Parameter()
    dataset = luigi.Parameter()
    meta_a = luigi.ListParameter()
    meta_b = luigi.ListParameter()
    meta_c = luigi.ListParameter()
    meta_x = luigi.ListParameter()
    blgsize = luigi.Parameter()

    def requires(self):
        return ConvertfBedToEigenstrat(self.group, self.dataset, GROUP_BY_POPS)

    def output(self):

        # shorten the meta pop names
        a = '-'.join([re.sub('[^A-Z]', '', meta) for meta in self.meta_a])
        b = '-'.join([re.sub('[^A-Z]', '', meta) for meta in self.meta_b])
        c = '-'.join([re.sub('[^A-Z]', '', meta) for meta in self.meta_c])
        x = '-'.join([re.sub('[^A-Z]', '', meta) for meta in self.meta_x])

        return [luigi.LocalTarget("qpf4ratio/{0}.{1}.a-{2}.b-{3}.c-{4}.x-{5}.blgsize-{6}.{7}".format(
            self.group, self.dataset, a, b, c, x, self.blgsize, ext)) for ext in ['par', 'log', 'poplist']]

    def run(self):

        # resolve the meta pops into lists of actual populations
        a_pops = get_metapops(self.group, self.dataset, self.meta_a)
        b_pops = get_metapops(self.group, self.dataset, self.meta_b)
        c_pops = get_metapops(self.group, self.dataset, self.meta_c)
        x_pops = get_metapops(self.group, self.dataset, self.meta_x)

        # get the outgroup
        o = OUTGROUP_POP[self.group] if self.group in OUTGROUP_POP else OUTGROUP_POP[self.dataset]

        # write the list of F4 ratio tests
        with self.output()[2].open('w') as fout:
            for a, b, c, x in itertools.product(a_pops, b_pops, c_pops, x_pops):
                # skip any duplicates
                if len(set([a, b, c, x])) == 4:
                    # f4(A,O; X,C) / f4(A,O; B,C)
                    fout.write("{a} {o} : {x} {c} :: {a} {o} : {b} {c}".format(a=a, b=b, c=c, x=x, o=o) + "\n")

        # compose the config settings
        config = [
            "genotypename: {}".format(self.input()[1].path),
            "snpname:      {}".format(self.input()[2].path),
            "indivname:    {}".format(self.input()[3].path),
            "popfilename:  {}".format(self.output()[2].path),
            "blgsize:      {}".format(self.blgsize)
        ]

        # the params need to be defined in a .par file
        parfile = self.output()[0].path

        with open(parfile, 'w') as par:
            par.write("\n".join(config))

        # run qp3pop
        log = run_cmd(["qpF4ratio", "-p", parfile])

        # save the log file
        with self.output()[1].open('w') as logfile:
            logfile.write(log)


class CTVTqpGraphPipeline(luigi.WrapperTask):
    """
    Run these specific qpGraph tasks from the CTVT pipeline
    """

    def requires(self):

        # new analysis group for Laurent mirroring the treemix group
        # yield QPGraphCluster('graph-pops2', 'merged_v2_TV_laurent', exhaustive=True)
        yield QPGraphCluster('graph-pops3', 'merged_v3_TV_laurent', exhaustive=True)


class CTVTFiguresPipeline(luigi.WrapperTask):
    """
    Run the specific elements of the CTVT pipeline
    """

    def requires(self):

        # Figure_NJTREE     / all-pops.merged_v3.njtree.pdf
        yield NeighborJoiningTree('all-pops', 'merged_v3')

        # Figure_NJVIET     / nj-pops.merged_v3_njviet.njtree.pdf
        yield NeighborJoiningTree('nj-pops', 'merged_v3_njviet')

        # Figure_PCA1       / all-pops.merged_v3.prj-DPC.PCA.1.2.pdf
        yield SmartPCAPlot('all-pops', 'merged_v3', ['DPC'], [(1, 2)])

        # Figure_PCA2       / dog-ctvt.merged_v3.prj-DPC.PCA.1.2.pdf
        yield SmartPCAPlot('dog-ctvt', 'merged_v3', ['DPC'], [(1, 2)])

        # Figure_PCA3       / dog-ctvt.merged_v3.prj-DPC-CTVT.PCA.1.2.pdf
        yield SmartPCAPlot('dog-ctvt', 'merged_v3', ['DPC', 'CTVT'], [(1, 2)])

        # Figure_TREEMIX    / graph-pops3.merged_v3_TV_laurent.treemix.geno.grp-pops.m0.pdf
        # Figure_TREEMIX1   / graph-pops3.merged_v3_TV_laurent.treemix.geno.grp-pops.m1.pdf
        # Figure_TREEMIX2   / graph-pops3.merged_v3_TV_laurent.treemix.geno.grp-pops.m2.pdf
        for m in range(0, TREEMIX_MAX_M + 1):
            # yield TreemixPlotM('graph-pops2', 'merged_v2_TV_laurent', GROUP_BY_POPS, m)
            yield TreemixPlotM('graph-pops3', 'merged_v3_TV_laurent', GROUP_BY_POPS, m)


class CTVTStatsPipeline(luigi.WrapperTask):
    """
    Run these specific qpGraph tasks from the CTVT pipeline
    """

    def requires(self):

        # SmartPCA
        for dataset in ['merged_v3', 'merged_v3_TV']:
            yield SmartPCAPlot('all-pops', dataset)
            yield SmartPCAPlot('dog-ctvt', dataset)
            yield SmartPCAPlot('dog-ctvt', dataset, ['DPC', 'CTVT'])

        yield SmartPCAPlot('all-pops', 'merged_SNParray_v5')

        # QPDstat
        for dataset in ['merged_v3', 'merged_v3_TV']:
            for blgsize in [1, 2]:
                yield QPDstat('all-pops', dataset, blgsize)

        # QPF4ratio
        a = ['CTVT']
        b = ['Pre-contact Dogs']
        c = ['European Dogs', 'East Asian Dogs']
        x = ['American Dogs']

        for dataset in ['merged_v3_hq', 'merged_v3_TV_hq', 'merged_SNParray_v5']:
            for blgsize in [1, 2]:
                yield QPF4ratio('all-pops', dataset, a, b, c, x, blgsize)

        # QP3Pop
        for dataset in ['merged_v3_hq', 'merged_v3_TV_hq']:
            yield QP3Pop('all-pops', dataset)


if __name__ == '__main__':
    luigi.run()
