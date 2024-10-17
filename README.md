# Pipeline for processing 1000 Genomes Phase 3 data sets

This is a `snakemake` pipeline for downloading and processing the [1000 Genomes Project (1kGP) Phase 3 data sets](https://www.internationalgenome.org/category/phase-3/) for use in genomic analyses. Ideally, you would integrate this into another workflow with the use of `snakemake`'s `module` statement, although you can use it in its own right, too. The workflow contains rules which will:
* download the raw sequence data and verify their checksums
* filter samples by ancestry
* remove related samples
* filter variants to retain SNPs only
* filter variants by MAF
* filter variants by QC parameters: genotype missingness, sample missingness, HWE test statistic
* filter variants to remove AT/GC SNPs
* prune variants by LD
* remove the MHC

Of course, there's much more you can do with these data and `plink`, but I intend this workflow to be a 'stub' of sorts which can be extended in the downstream direction with additional rules after import as a `module`. Hopefully these rules suffice to get you started.

The `assembly` wildcard allows you to specify either the `hg19` or `hg38` reference genome. You should prefer the latter where possible as the `hg38` data were generated through [high-coverage resequencing of the 1kGP samples](https://doi.org/10.1016/j.cell.2022.08.004), and the `hg38` is itself a more accurate assembly and seven years old at this point; `hg19` is ancient, stop using it if you can! Check out `workflow/rules/wildcard_constraints.smk` for other wildcards, many of which can only take a restricted set of values (like an enum).

Dependencies for the software are provided in a `docker` container [hosted on DockerHub](https://hub.docker.com/repository/docker/twillis209/1kgp-pipeline/general). `snakemake` should pull this automatically so there's no need to build the container yourself, although I have included the Dockerfile under the `docker` directory so you can see what it contains. The image is derived from the [`rocker/r-ver` image](https://rocker-project.org/images/versioned/r-ver).


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
This allows for more parsimonious rules in the `smk` files. In my experience, the `plink`-based rules will generally run faster the more you pump up the thread count.

# Outstanding issues (17/10/24)

I am able to run this locally, but the cluster I use seems not to allow `plink` to create log files. I hope this is just some idiosyncrasy of cluster configuration, but it may relate to `apptainer`'s interaction with the file system, something that can be configured via the `snakemake` CLI. I'm in the process of troubleshooting this.

If you happen upon this repo and have a problem, please open an issue.
