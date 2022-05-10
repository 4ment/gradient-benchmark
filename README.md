# autodiff-experiments

[![Docker Image CI](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml/badge.svg)](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml)

This benchmark compares the efficiency and accuracy of ...


| Program      | Language  | Framework    | Gradient | [BITO] |
| ------------ | --------- | ------------ | ---------| -----|
| [physher]    | C         |              | analytic |     |
| [phylostan]  | Stan      | [Stan]       | autodiff |     |
| [phylojax]   | python    | [JAX]        | autodiff |     |
| [torchtree] | python    | [PyTorch]    | autodiff | yes |
| [treeflow]   | python    | [TensorFlow] | autodiff | yes |

The gradient of the tree likelihood is optionaly computed by [BITO], an efficient C++ library that analytically calculate
the gradient using the [BEAGLE] library. BITO is only available in treeflow and torchtree.

## Dependencies
You will need to install [nextflow](https://www.nextflow.io) and [docker](https://www.docker.com) to run this benchmark.
Docker is not required but it is highly recommended to use it due to the numerous dependencies.

## Installation

    git clone 4ment/autodiff-experiments.git

### Initialize treetime_validation

    git submodule update --init --recursive

## Running the pipeline with docker

    nextflow run 4ment/autodiff-experiments -profile docker

[physher]: https://github.com/4ment/physher
[phylostan]: https://github.com/4ment/phylostan
[phylojax]: https://github.com/4ment/phylojax
[torchtree]: https://github.com/4ment/torchtree
[treeflow]: https://github.com/christiaanjs/treeflow

[BITO]: https://github.com/phylovi/bito
[BEAGLE]: https://github.com/beagle-dev/beagle-lib

[Stan]: https://mc-stan.org
[JAX]: https://github.com/google/jax
[PyTorch]: https://pytorch.org
[TensorFlow]: https://www.tensorflow.org