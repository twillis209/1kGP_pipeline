# Pipeline for processing 1000 Genomes Phase 3 data sets

This is a `snakemake` pipeline for downloading and processing the [1000 Genomes Project (1kGP) Phase 3 data sets](https://www.internationalgenome.org/category/phase-3/) for use in genomic analyses. Ideally, you would integrate this into another workflow with the use of `snakemake`'s `module` statement, although you can use it in its own right, too. The workflow contains rules which will:
* download the raw sequence data in `vcf.gz` format and verify their checksums
* filter samples by ancestry
* remove related samples
* filter variants to retain SNPs only
* filter variants by MAF
* filter variants by QC parameters: genotype missingness, sample missingness, HWE test statistic
* filter variants to remove AT/GC SNPs
* prune variants by LD
* remove the MHC
* add the centimorgan field to `plink`-style bfiles for the calculation of LD scores
* compute LD scores

Of course, there's much more you can do with these data and `plink`, but I intend this workflow to be a 'stub' of sorts which can be extended in the downstream direction with additional rules after import as a `module`. Hopefully these rules suffice to get you started.

If you're use the pipeline for the first time, I recommend testing it with chromosome 22 (`chr22`) which, as the smallest chromosome, should allow a quick test run.

The `assembly` wildcard allows you to specify either the `hg19` or `hg38` reference genome. You should prefer the latter where possible as the `hg38` data were generated through [high-coverage resequencing of the 1kGP samples](https://doi.org/10.1016/j.cell.2022.08.004), and the `hg38` is itself a more accurate assembly and seven years old at this point; `hg19` is ancient, stop using it if you can! Check out `workflow/rules/wildcard_constraints.smk` for other wildcards, many of which can only take a restricted set of values (like an enum) or must match a regular expression.

Dependencies for the software are provided in a `docker` container [hosted on DockerHub](https://hub.docker.com/repository/docker/twillis209/1kgp-pipeline/general). `snakemake` should pull this automatically so there's no need to build the container yourself, although I have included the Dockerfile under the `docker` directory so you can see what it contains. The image is derived from the [`rocker/r-ver` image](https://rocker-project.org/images/versioned/r-ver).

`pandas` is required for some `run` rules and this is provided in the `conda` environment defined in `envs/global.yaml`. Your profile should contain the following so that both `docker` and `conda` are used:
```
software-deployment-method: "apptainer"
software-deployment-method: "conda"
```
At the moment, the `conda` statement required to inject this dependency into the `snakemake` process running the workflow is omitted as it appears not to play nicely with the `module` statement I use to import this workflow into other workflows (e.g. the [`GWAS_tools`](https://github.com/twillis209/GWAS_tools)), so I'm afraid you'll have to have `pandas` in whichever environment `snakemake` is running in.

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

## Importing this workflow in another `snakemake` workflow using `module`

I wrote this workflow so it could be plugged into others with use of the `module` statement. I do this like so in the [`GWAS_tools` worfklow](https://github.com/twillis209/GWAS_tools):

```
module kGP_pipeline:
    snakefile: github("twillis209/1kGP_pipeline", path = "workflow/Snakefile", commit = 'master')
    config: config['1kGP_pipeline']
```

The config handling is currently a bit hacky: rather than the module reading its own config file (i.e. the one in `1kGP_pipeline/workflow/config/config.yaml`), I have to copy the contents of that file to the config file for `GWAS_tools` and reference that when overwriting the `1kGP_pipeline` config with the `config:` statement. This doesn't scale very well when importing `GWAS_tools` into still larger workflows as I like to do, but I'm optimistic I can find a better solution within `snakemake`.

## A version without containers

At present I'm not able to resolve the issue I describe below relating to log files in my current HPC environment, so have had to roll back the containerisation in the dedicated `conda` branch.

## LD score calculation

Including this functionality in the pipeline is arguably scope creep of the type I'd like to avoid, but LD scores are near-ubiquitous these days across analyses downstream of GWAS. Note that I've followed the guidance on the [LDSC GitHub](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) when calculating the scores, but a more careful approach would probably adjust the window parameter for the MHC and any other regions featuring long-range linkage disequilibrium patterns. I can't guarantee that is sufficient for dealing with the MHC, though: for immune-related phenotypes it contains variants with extremely large effects which may not be handled well by the heritability estimation methods like LDSC or SumHer (the authors of the latter suggest inclusion of fixed effects to account for this).

`plink` does not like inclusion of the pseudoautosomal regions on the X chromosome (it won't work with chromosomes labelled 'PAR1' and 'PAR2'), so I drop these for calculation of LD scores (see the `sans_pars` value for the `variant_set` wildcard.

# Outstanding issues (10/2/25)

I am able to run this locally, but the cluster I use seems not to allow `plink` to create log files. I hope this is just some idiosyncrasy of cluster configuration, but it may relate to `apptainer`'s interaction with the file system, something that can be configured via the `snakemake` CLI. I'm in the process of troubleshooting this and apparently [I'm not the only one with this issue](https://github.com/snakemake/snakemake/issues/2959).

The `docker` branch is lagging behind quite a bit; the `conda` branch is the current branch and the only one to include the LD score functionality.

If you happen upon this repo and have a problem, please open an issue.
