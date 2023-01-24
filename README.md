# autodiff-experiments

[![Docker Image CI](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml/badge.svg)](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml)

This repository contains the pipeline and data sets supporting the results of the following article:

Fourment M, Swanepoel CJ, Galloway JG, Ji X, Gangavarapu K, Suchard MA, Matsen IV FA. Automatic differentiation is no panacea for phylogenetic gradient computation. [arXiv:2211.02168](https://arxiv.org/abs/2211.02168)

This benchmark compares the efficiency (memory usage and speed) of several gradient implementations of phylogenetic models (e.g., tree likelikelihood and coalescent model). The goal of this study is to compare the efficiency of automatic differentiation (AD) and analytic gradient. The pipeline reuses parts of the treetime validation [workflow](https://github.com/neherlab/treetime_validation).


| Program      | Language  | Framework    | Gradient | [BITO] support |
| ------------ | --------- | ------------ | :-------:| :-----:|
| [physher]    | C         |              | analytic |     |
| [phylostan]  | Stan      | [Stan]       | AD       |     |
| [phylojax]   | python    | [JAX]        | AD       |     |
| [torchtree]  | python    | [PyTorch]    | AD       | :white_check_mark: |
| [treeflow]   | python    | [TensorFlow] | AD       |  |

The gradient of the tree likelihood is optionaly computed by [BITO], an efficient C++ library that analytically calculate
the gradient using the [BEAGLE] library. torchtree uses the [torchtree-bito] plugin to access BITO.

## Dependencies
You will need to install [nextflow](https://www.nextflow.io) and [docker](https://www.docker.com) to run this benchmark.
Docker is not required but it is highly recommended to use it due to the numerous dependencies.

## Installation

    git clone 4ment/autodiff-experiments.git

### Initialize treetime_validation

    git submodule update --init --recursive

## Running the pipeline with docker

    nextflow run 4ment/autodiff-experiments -profile docker -with-trace

Since the pipeline will take weeks to run to completion one should use a high performance computer. Examples of configuration files for pbspro and slurm can be found in the [configs](configs/) folder.

## Summarizing results

Before generating the figures, we need to extract memory usage information from the `trace.txt` file and `work` directory:

    python scripts/parse-trace.py work/ trace.txt > results/trace.csv

Generate figures in a single pdf:

    Rscript -e 'rmarkdown::render("plot.Rmd")'

## Library versions

For reproducbility, we provide below the version or commit hash of each library/program used in the benchmark.

| Library      | Version  |
| ------------ | -------- |
| [jax]                  | 0.2.24   |
| jaxlib                 | 0.3.7    |
| numpy                  | 1.22     |
| [pystan]               | 2.19.1.1 |
| [tensorflow]           | 2.10.0   |
| tensorflow-probability | 0.18.0   |
| [pytorch]              | 1.12.1   |

 | Program | Version/hash |
 | ------------ | -------- |
 | [bito] |         cc0806abcd0b9f2fab604e800c674c9a5c5afebe |
 | [phylojax]       | a1612cae36292af76e8d24cc40d6544162c987aa |
 | [phylostan]      | 1.0.5 |
 | [physher]        | b19ff2f9422f29ba1ab31306a3fe29ab6a6f607b |
 | [torchtree]      | f3831650a807e74cc2e9478009e57a41f47bed8d |
 | [torchtree-bito] | e2a95cefb13968f95f6e5520bd0a52d726ee7fc9 |
 | [treeflow]       | e3414dcc9e764d06abc3e19c1d0f55110499e2ea |
 


[physher]: https://github.com/4ment/physher
[phylostan]: https://github.com/4ment/phylostan
[phylojax]: https://github.com/4ment/phylojax
[torchtree]: https://github.com/4ment/torchtree
[torchtree-bito]: https://github.com/4ment/torchtree-bito
[treeflow]: https://github.com/christiaanjs/treeflow

[BITO]: https://github.com/phylovi/bito
[BEAGLE]: https://github.com/beagle-dev/beagle-lib

[Stan]: https://mc-stan.org
[JAX]: https://github.com/google/jax
[PyTorch]: https://pytorch.org
[TensorFlow]: https://www.tensorflow.org
[pystan]: https://github.com/stan-dev/pystan2