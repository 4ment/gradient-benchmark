# autodiff-experiments

[![Docker Image CI](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml/badge.svg)](https://github.com/4ment/autodiff-experiments/actions/workflows/docker-image.yml)

This benchmark compares the efficiency and accuracy of ...


| Program      | Language  | Framework |
| ------------ | --------- | --------- |
| [physher]    | C         | |
| [phylostan]  | Stan      | [Stan](https://mc-stan.org) |
| [phylojax]   | python    | [JAX](https://github.com/google/jax) |
| [phylotorch] | python    | [PyTorch](https://pytorch.org) |
| [treeflow]   | python    | [TensorFlow](https://www.tensorflow.org) |

## Dependencies
You will need to install [nextflow](https://www.nextflow.io) and [docker](https://www.docker.com) to run this benchmark.
Docker is not required but it is highly recommended to use it due to the numerous dependencies.

## Running the pipeline with docker

    nextflow run 4ment/autodiff-experiments -profile docker

[physher]: https://github.com/4ment/physher
[phylostan]: https://github.com/4ment/phylostan
[phylojax]: https://github.com/4ment/phylojax
[phylotorch]: https://github.com/4ment/phylotorch
[treeflow]: https://github.com/christiaanjs/treeflow