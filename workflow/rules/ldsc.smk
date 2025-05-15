rule write_per_chrom_bfiles_for_ld_score_estimation:
    input:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm", ".bed", ".bim", ".fam")
    output:
        temp(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm/chr{chr_no}", ".bed", ".bim", ".fam"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.bed')
    threads: 8
    conda: env_path("global.yaml")
    shell: "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --chr {wildcards.chr_no} --make-bed --out {params.out_stem}"

rule compute_ld_scores:
    input:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm/chr{chr_no}", ".bed", ".bim", ".fam")
    output:
        multiext("results/ldsc/ld_scores/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/chr{chr_no}", ".l2.M", ".l2.ldscore.gz", ".l2.M_5_50")
    log:
        "results/ldsc/ld_scores/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/chr{chr_no}.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.l2.M')
    threads: 8
    resources:
        runtime = 10,
    conda: env_path("ldsc.yaml")
    shell: "ldsc.py --bfile {params.in_stem} --l2 --ld-wind-cm 1 --out {params.out_stem}"

rule compute_all_ld_scores:
    input:
       [f"results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr{x}.l2.M" for x in range(1,23)]
