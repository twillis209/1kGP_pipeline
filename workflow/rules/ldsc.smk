rule write_out_per_chrom_recombination_map_files:
    input:
        "resources/1kG/{assembly}/genetic_map_{assembly}_withX.txt.gz"
    output:
        temp([f"resources/1kG/{{assembly}}/genetic_map_{{assembly}}/chr{x}.txt" for x in list(range(1,23))+['X']])
    params:
        root_dir = subpath(input[0], strip_suffix = "_withX.txt.gz")
    localrule: True
    threads: 1
    resources:
        runtime = 20
    shell:
        """
        for x in {{1..22}}; do
            echo -e "position\trrate\tgposition" >"{params.root_dir}/chr"$x.txt
        done

        echo -e "position\trrate\tgposition" >"{params.root_dir}/chrX.txt"

        zcat {input} | tail -n +2 | awk 'BEGIN {{OFS="\t"}} {{print $2, $3, $4 >> "{params.root_dir}/chr"$1".txt"}}'
        """

rule write_out_bed_format_files_with_cm_field:
    input:
        rules.write_out_qced_data_to_bed_format.output,
        map_files = [f"resources/1kG/{{assembly}}/genetic_map_{{assembly}}/chr{x}.txt" for x in list(range(1,23))+['X']]
    output:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm", ".bed", ".bim", ".fam")
    log:
        log_file = "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm.log"
    params:
        in_stem = subpath(input[0], strip_suffix = ".bed"),
        out_stem = subpath(output[0], strip_suffix = '.bed'),
        map_pattern = "resources/1kG/{assembly}/genetic_map_{assembly}/chr@.txt"
    threads: 16
    resources:
        runtime = 5
    group: "1kG"
    conda: env_path("global.yaml")
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --cm-map {params.map_pattern} --make-bed --out {params.out_stem} >{log.log_file}"

rule write_per_chrom_bfiles_for_ld_score_estimation:
    input:
        rules.write_out_bed_format_files_with_cm_field.output
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
        rules.write_per_chrom_bfiles_for_ld_score_estimation.output
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
