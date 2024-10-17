# Pipeline for processing 1000 Genomes Phase 3 data sets

This is a `snakemake` pipeline for downloading and processing the [1000 Genomes Project (1kGP) Phase 3 data sets](https://www.internationalgenome.org/category/phase-3/) for use in genomic analyses. Ideally, you would integrate this into another workflow with use of `snakemake`'s `module` statement, although you can use it in its own right, too. You'll need a 

Dependencies for the software are provided in a `docker` container [hosted on DockerHub](https://hub.docker.com/repository/docker/twillis209/1kgp-pipeline/general). `snakemake` should pull this automatically so there's no need to build the container yourself, although I have included the Dockerfile under the `docker` directory if you're interested in what it contains.

TODO:
* manage `pandas` and other dependencies in `envs/global.yaml`, if `run` is used in a rule the container is not invoked (AFAIK)
