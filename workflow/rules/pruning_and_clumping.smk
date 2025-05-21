rule create_pruned_ranges:
    input:
        rules.filter_variant_set.output,
    output:
        temp(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out"))
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.prune.in'),
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    resources:
        runtime = 20
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --indep-pairwise {wildcards.window_size}kb 1 {params.r2} --out {params.out_stem}"

rule prune_variants:
    input:
        rules.filter_variant_set.output,
        range_file = rules.create_pruned_ranges.output[1]
    output:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged", ".pgen", ".psam", ".pvar.zst")
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen')
    threads: 16
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --exclude {input.range_file} --make-pgen vzs --out {params.out_stem}"

# NB: plink reports the frequency of the alternative allele, not the MAF
rule compute_alt_freqs:
    input:
        rules.filter_variant_set.output,
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged.afreq"
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.afreq')
    threads: 16
    resources:
        runtime = 10
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --freq --out {params.out_stem}"

rule compute_maf_clumping_statistic:
    """
    NB: the smaller the statistic, the more plink prefers the SNP when clumping, so we use major, rather than minor, allele frequency, to clump with a preference for higher minor allele frequency SNPs
    """
    input:
        rules.compute_alt_freqs.output
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/clumping_stat.tsv"
    run:
        daf = pl.read_csv(input[0], separator = '\t', schema_overrides = {"#CHROM": pl.String()})

        daf = daf.with_columns(
            pl.when(pl.col("ALT_FREQS") > 0.5)
            .then(1 - pl.col("ALT_FREQS"))
            .otherwise(pl.col("ALT_FREQS"))
            .alias("MINOR_AF"),
        )

        daf = daf.with_columns(
            (1 - pl.col("MINOR_AF")).alias("MAJOR_AF")
        )

        daf.write_csv(output[0], separator = '\t')

rule clump_by_maf:
    input:
        rules.filter_variant_set.output,
        stat = rules.compute_maf_clumping_statistic.output
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/clumped_by_maf/{window_size}_{r2}/merged.clumps"
    log:
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen'),
        clump_field = "MAJOR_AF",
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    group: "1kG"
    shell: "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs  --out {params.out_stem} --clump-field {params.clump_field} --clump {input.stat} --clump-r2 {params.r2} --clump-kb {wildcards.window_size} --clump-p1 1"
