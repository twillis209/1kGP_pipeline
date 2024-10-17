# Pipeline for processing 1000 Genomes Phase 3 data sets

This is a `snakemake` pipeline for downloading and processing the [1000 Genomes Project (1kGP) Phase 3 data sets](https://www.internationalgenome.org/category/phase-3/) for use in genomic analyses. Ideally, you would integrate this into another workflow with use of `snakemake`'s `module` statement, although you can use it in its own right, too. The `assembly` wildcard allows you to specify the `hg19` or `hg38` reference genomes. You should prefer the latter where possible as the `hg38` data were generated through [high-coverage resequencing of the 1kGP samples](https://doi.org/10.1016/j.cell.2022.08.004), and the `hg38` is itself a more accurate assembly and seven years old at this point; `hg19` is ancient, stop using it!

Dependencies for the software are provided in a `docker` container [hosted on DockerHub](https://hub.docker.com/repository/docker/twillis209/1kgp-pipeline/general). `snakemake` should pull this automatically so there's no need to build the container yourself, although I have included the Dockerfile under the `docker` directory if you're interested in what it contains.

`pandas` is required for some `run` rules and this is provided in the `conda` environment defined in `envs/global.yaml`. Your profile should contain the following so that both `docker` and `conda` are used:
```
software-deployment-method: "apptainer"
software-deployment-method: "conda"
```

I don't define the `mem_mb` resource for each rule as that number is highly platform-specific, but I do recommend setting it as follows in your profile:
```
default-resources:
  mem_mb: threads * <RAM per core>
```
For example, the Intel 'Cascade Lake' CPUs I use on my local cluster provide 3420MB per core, so I set mine to:
```
default-resources:
  mem_mb: threads * 3420
```
This gives you much cleaner, more parsimonious rules.

In my experience, the `plink`-based rules will generally run faster the more you pump up the thread count.
