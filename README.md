# CTVT
This repository contains Python and R code for the ancestry analyses of the ancient nuclear DNA from the paper 
[The evolutionary history of dogs in the Americas](https://doi.org/10.1126/science.aao4776).

![Figure qpGraph](./Figure_QPGRAPH.png?raw=true)

## Citation
If you reuse any of this code then please cite the paper:
> Leathlobhair, M.N., Perri, A.R., Irving-Pease, E.K., Witt, K.E., Linderholm, A., Haile, J., Lebrasseur, O., Ameen, C., 
> Blick, J., Boyko, A.R., Brace, S., Cortes, Y.N., Crockford, S.J., Devault, A., Dimopoulos, E.A., Eldridge, M., Enk, 
> J., Gopalakrishnan, S., Gori, K., Grimes, V., Guiry, E., Hansen, A.J., Hulme-Beaman, A., Johnson, J., Kitchen, A., 
> Kasparov, A.K., Kwon, Y.-M., Nikolskiy, P.A., Lope, C.P., Manin, A., Martin, T., Meyer, M., Myers, K.N., Omura, M., 
> Rouillard, J.-M., Pavlova, E.Y., Sciulli, P., Sinding, M.-H.S., Strakova, A., Ivanova, V.V., Widga, C., Willerslev, 
> E., Pitulko, V.V., Barnes, I., Gilbert, M.T.P., Dobney, K.M., Malhi, R.S., Murchison, E.P., Larson, G., Frantz, 
> L.A.F., 2018. The evolutionary history of dogs in the Americas. *Science* 361, 81–85. 
> https://doi.org/10.1126/science.aao4776

## Installation

To reproduce the analyses from the paper you will need to install the following dependencies.

### Python

Python ≥ 2.7 with the following modules:

* [dendropy](https://github.com/jeetsukumaran/DendroPy)
* [graphviz](https://github.com/xflr6/graphviz)
* [luigi](https://github.com/spotify/luigi)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [numpy](https://github.com/numpy/numpy)
* [pathos](https://github.com/uqfoundation/pathos)
* [psutil](https://github.com/giampaolo/psutil)
* [pyparsing](https://github.com/pyparsing/pyparsing)
* [scipy](https://github.com/scipy/scipy)


```bash
pip install dendropy graphviz luigi matplotlib numpy pathos psutil pyparsing scipy 
```

The full list of Python modules installed in the project environment can be
found in the [requirements](requirements.txt) file.

### R

R ≥ 3.4 with the following modules:

* [ape](https://cran.r-project.org/web/packages/ape/)
* [combinat](https://cran.r-project.org/web/packages/combinat/)
* [doParallel](https://cran.r-project.org/web/packages/doParallel)
* [foreach](https://cran.r-project.org/web/packages/foreach)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/)
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/)
* [scales](https://cran.r-project.org/web/packages/scales/)
* [stringr](https://cran.r-project.org/web/packages/stringr/)


```R
install.packages(c("ape", "combinat", "doParallel", "foreach", "ggplot2", "gridExtra", "RColorBrewer", "reshape2", 
                   "scales", "stringr"))
```

### Other

* [AdmixTools](https://github.com/DReichLab/AdmixTools)
* [ADMIXTURE](http://software.genetics.ucla.edu/admixture/)
* [graph_tool](https://graph-tool.skewed.de/)
* [Plink](http://www.cog-genomics.org/plink2)
* [TreeMix](https://bitbucket.org/nygcresearch/treemix/)

## Running the pipeline

The pipeline is broken into three [luigi](https://github.com/spotify/luigi) wrapper tasks.

### Stats pipeline

* run qpDstat
* run qpF4ratio
* run qp3Pop

```bash
luigi --module pipeline_ctvt CTVTStatsPipeline
```

### qpGraph pipeline
 
* fit qpGraph models

```bash
luigi --module pipeline_ctvt CTVTqpGraphPipeline
```

### Supplementary figures pipeline

* run SmartPCA
* plot NJ trees
* run Treemix


```bash
luigi --module pipeline_ctvt CTVTFiguresPipeline
```


## Author

Evan K. Irving-Pease, [PalaeoBARN](https://www.palaeobarn.com/), University of Oxford 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
