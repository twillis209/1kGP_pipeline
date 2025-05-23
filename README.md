# Pipeline for processing the 1000 Genomes Phase 3 datasets

This is a `snakemake` pipeline for downloading and processing the [1000 Genomes Project (1kGP) Phase 3 datasets](https://www.internationalgenome.org/category/phase-3/) for use in genomic analyses. Ideally, you would integrate this into another workflow with the use of `snakemake`'s `module` statement, although you can use it in its own right, too. The workflow contains rules which will:
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

Dependencies like `plink`, `plink2`, and `pandas` are managed with `conda`. `snakemake` should build the environments automatically, although in my experience it's quite slow in doing so compared with running `conda` in other contexts.

An earlier version of the pipeline used `docker` instead, but I no longer maintain the `docker` branch on which this resides. The image was derived from the [`rocker/r-ver` image](https://rocker-project.org/images/versioned/r-ver).

`pandas` is required for some `run` rules and this is provided in the `conda` environment defined in `envs/global.yaml`.

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

The config handling is currently a bit hacky: rather than the module reading its own config file (i.e. the one in `1kGP_pipeline/workflow/config/config.yaml`), I have to copy the contents of that file to the config file for `GWAS_tools` and reference that when overwriting the `1kGP_pipeline` config with the `config:` statement. This doesn't scale very well when importing `GWAS_tools` into still larger workflows as I like to do.

Re: `GWAS_tools` in its own right, I don't use or recommend it anymore and now use the EBI's [EBISPOT harmoniser](https://github.com/EBISPOT/gwas-sumstats-harmoniser).

## LD score calculation

Including this functionality in the pipeline is arguably scope creep of the type I'd like to avoid, but LD scores are near-ubiquitous these days across analyses downstream of GWAS. Note that I've followed the guidance on the [LDSC GitHub](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) when calculating the scores, but a more careful approach would probably adjust the window parameter for the MHC and any other regions featuring long-range linkage disequilibrium patterns. I can't guarantee that is sufficient for dealing with the MHC, though: for immune-related phenotypes it contains variants with extremely large effects which may not be handled well by the heritability estimation methods like LDSC or SumHer (the authors of the latter suggest inclusion of fixed effects to account for this). You can include the token `sans_long_range_ld` in the `variant_set` wildcard to have `plink` remove such regions (see below for more on these regions).

`plink` does not like inclusion of the pseudoautosomal regions on the X chromosome (it won't work with chromosomes labelled 'PAR1' and 'PAR2'), so I drop these for calculation of LD scores: inclusion of the `sans_pars` token in the `variant_set` wildcard should work.

# Acknowledgements

## Long-range LD regions

The long-range LD region bedfiles are taken from the Meyer lab's [`plinkQC`](https://github.com/meyer-lab-cshl/plinkQC) package. There's a page on these regions, including a nice depiction of their positions in the genome, at the [Abecasis group's wiki](https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)). That page cites Carl Anderson's [2010 article](https://www.nature.com/articles/nprot.2010.116) for these regions, but the article is paywalled and the PubMed version doesn't seem to contain a table of the regions.
