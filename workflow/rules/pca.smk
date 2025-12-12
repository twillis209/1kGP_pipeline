rule pca:
    input:
        rules.prune_variants.output
    output:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/{window_size}_1_{r2}/pca/merged", eigenvec = ".eigenvec", eigenval = ".eigenval", allele = ".eigenvec.allele", acount = ".acount")
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output.eigenvec, strip_suffix = '.eigenvec'),
        no_of_pcs = 20
    threads: 32
    shell: """
        plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --freq counts --pca {params.no_of_pcs} allele-wts vcols=chrom,ref,alt --out {params.out_stem}
        """

rule score_samples_with_pca:
    """
    This uses the newer plink2 `score` program, see https://www.cog-genomics.org/plink/2.0/score#pca_project. In particular note the following from the docs: '[T]hese PCs will be scaled a bit differently from ref_data.eigenvec; you need to multiply or divide the PCs by a multiple of sqrt(eigenvalue) to put them on the same scale.'

    """
    input:
        rules.pca.output,
        rules.prune_variants.output
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/{window_size}_1_{r2}/pca/merged.sscore"
    params:
        pfile_stem = subpath(input[4], strip_suffix = '.pgen'),
        pca_stem = subpath(input[0], strip_suffix = '.eigenvec'),
        out_stem = subpath(output[0], strip_suffix = '.sscore'),
        variant_id_col = 2,
        allele_code_col = 5,
        score_col_range = "6-25",
    threads: 16
    shell:
        """
        plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.pfile_stem} vzs --read-freq {params.pca_stem}.acount --score {params.pca_stem}.eigenvec.allele {params.variant_id_col} {params.allele_code_col}  header-read no-mean-imputation variance-standardize --score-col-nums {params.score_col_range} --out {params.out_stem}
        """

rule plot_scores_on_first_two_pcs:
    input:
        scores = rules.score_samples_with_pca.output,
        ped = "resources/1kG/hg38/ped.txt"
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/{window_size}_1_{r2}/pca/merged.sscore.png"
    localrule: True
    conda: "../envs/r.yaml"
    script: "../scripts/plot_scores_on_first_two_pcs.R"
