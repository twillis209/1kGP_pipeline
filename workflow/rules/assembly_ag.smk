rule make_1kG_sex_file:
    input:
        "resources/1kG/{assembly}/ped.txt",
    output:
        "results/1kG/{assembly}/sex.tsv"
    localrule: True
    run:
        ped = pd.read_csv(input[0], sep = '\t', header = 0)

        if wildcards.assembly == 'hg19':
            ped = ped[['Family ID', 'Individual ID', 'Gender']]
            ped = ped.rename({'Family ID': '#FID', 'Individual ID': 'IID', 'Gender': 'Sex'}, axis = 1)
        else:
            ped = ped[['FamilyID', 'SampleID', 'Sex']]
            ped = ped.rename({'FamilyID': '#FID', 'SampleID': 'IID'}, axis = 1)

        ped.to_csv(output[0], index = False, sep = '\t')

rule vcf_to_pgen:
    input:
        vcf = "resources/1kG/{assembly}/{chr}.vcf.gz",
        sex = "results/1kG/{assembly}/sex.tsv",
        ref = "resources/genome_reference/{assembly}.fa.zst"
    output:
        temp(multiext("results/1kG/{assembly}/{chr}", ".pgen", ".pvar.zst", ".psam"))
    params:
        out = "results/1kG/{assembly}/{chr}",
        id_format = "@:#:\$r:\$a",
        max_allele_len = 20
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input.vcf} --make-pgen vzs --fa {input.ref} --ref-from-fa 'force' --out {params.out} --set-all-var-ids {params.id_format} --max-alleles 2 --new-id-max-allele-len {params.max_allele_len} truncate --update-sex {input.sex} --split-par '{wildcards.assembly}'"

rule make_1kG_unrelated_sample_files:
     input:
         ped = "resources/1kG/{assembly}/ped.txt",
         # No hg38 panel file required
         panel = "resources/1kG/hg19/panel.txt",
     output:
         eur = "results/1kG/{assembly}/eur.samples",
         afr = "results/1kG/{assembly}/afr.samples",
         amr = "results/1kG/{assembly}/amr.samples",
         eas = "results/1kG/{assembly}/eas.samples",
         sas = "results/1kG/{assembly}/sas.samples",
         all = "results/1kG/{assembly}/all.samples"
     localrule: True
     conda:
         env_path("pid_cfdr_pipeline.yaml")
     script:
        script_path("write_1kG_sample_files.R")

rule get_ancestry_specific_samples:
     input:
        multiext("results/1kG/{assembly}/{chr}", ".pgen", ".pvar.zst", ".psam"),
        sample_file = "results/1kG/{assembly}/{ancestry}.samples"
     output:
        temp(multiext("results/1kG/{assembly}/{ancestry}/{chr}", ".pgen", ".pvar.zst", ".psam")),
     log:
        log = "results/1kG/{assembly}/{ancestry}/{chr}.log"
     params:
        in_stem = "results/1kG/{assembly}/{chr}",
        out_stem = "results/1kG/{assembly}/{ancestry}/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     group: "1kG"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --keep {input.sample_file} --make-pgen vzs --out {params.out_stem} &>{log.log}"

rule retain_snps_only:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/{chr}", ".pgen", ".pvar.zst", ".psam"),
    output:
        temp(multiext("results/1kG/{assembly}/{ancestry}/snps_only/{chr}", ".pgen", ".pvar.zst", ".psam"))
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/{chr}",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/{chr}",
        maf = 0.005
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --snps-only 'just-acgt' --rm-dup 'force-first' --maf {params.maf} --make-pgen vzs --out {params.out_stem}"

rule merge_pgen_files:
    input:
        expand("results/1kG/{{assembly}}/{{ancestry}}/snps_only/{chr}.{ext}", chr = [f"chr{x}" for x in range(1,23)]+["chrX"], ext = ["pgen", "pvar.zst", "psam"])
    output:
        protected(multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged", ".pgen", ".pvar.zst", ".psam")),
        pmerge_file = "results/1kG/{assembly}/{ancestry}/snps_only/pmerge.txt"
    params:
        in_dir = "results/1kG/{assembly}/{ancestry}/snps_only",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged",
        max_allele_ct = 2
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 10
    group: "1kG"
    shell: """
        for i in {{1..22}}; do
        echo "chr$i" >>{output.pmerge_file}
        done

        echo "chrX" >>{output.pmerge_file}

        plink2 --memory {resources.mem_mb} --threads {threads} --pmerge-list {output.pmerge_file} pfile-vzs --pmerge-list-dir {params.in_dir} --merge-max-allele-ct {params.max_allele_ct} --pmerge-output-vzs --out {params.out_stem}
    """

rule pgen_to_hap_and_legend:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/merged", ".haps", ".legend", ".sample"))
    params:
        stem = "results/1kG/{assembly}/{ancestry}/{variant_type}/merged"
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 120
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --export hapslegend --out {params.stem}"

rule compute_maf:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        "results/1kG/{assembly}/{ancestry}/{variant_type}/merged.afreq"
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/{variant_type}/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/{variant_type}/merged",
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 120
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --freq --out {params.out_stem}"

rule write_out_merged_bed_format_files:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged", ".pgen", ".pvar.zst", ".psam"),
    output:
        temp(multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged", ".bed", ".bim", ".fam"))
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged",
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 5
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --make-bed --out {params.in_stem}"

rule qc:
     input:
         multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged", ".pgen", ".pvar.zst", ".psam")
     output:
         multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/merged", ".pgen", ".pvar.zst", ".psam")
     log:
         log = "results/1kG/{assembly}/{ancestry}/snps_only/qc/merged.log"
     params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/merged",
        geno = 0.01,
        mind = 0.01,
        hwe = 1e-50
     threads: 16
     resources:
        mem_mb = get_mem_mb,
        runtime = 10
     group: "1kG"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --geno {params.geno} --mind {params.mind} --hwe {params.hwe} --make-pgen vzs --out {params.out_stem} &>{log.log}"

rule decompress_pvar_for_at_gc_snps:
    input:
        "results/1kG/{assembly}/{ancestry}/snps_only/qc/merged.pvar.zst"
    output:
        temp("results/1kG/{assembly}/{ancestry}/snps_only/qc/merged.pvar")
    threads: 1
    resources:
        runtime = 120
    shell:
        "zstdcat {input} | grep -v '^#' >{output}"

rule identify_at_gc_snps:
    input:
        "results/1kG/{assembly}/{ancestry}/snps_only/qc/merged.pvar"
    output:
        "results/1kG/{assembly}/{ancestry}/snps_only/qc/at_gc_snps.txt"
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = 60
    shell: """
    """

rule remove_at_gc_snps:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/merged", ".pgen", ".pvar.zst", ".psam"),
        at_gc_variants = "",
    output:
        multiext("results/1kG/{assembly}/{ancestry}/sans_at_gc_snps_only/qc/merged", ".pgen", ".pvar.zst", ".psam")
    log:
        log = "results/1kG/{assembly}/{ancestry}/sans_at_gc_snps_only/qc/merged.log"
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/sans_at_gc_snps_only/qc/merged"
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 10
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --exclude {input.at_gc_variants} --make-pgen vzs --out {params.out_stem} &>{log.log}"

rule copy_to_all_variant_set:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/all/merged", ".pgen", ".pvar.zst", ".psam")
    params:
        out = "results/1kG/{assembly}/{ancestry}/snps_only/qc/all"
    threads: 1
    resources:
        mem_mb = get_mem_mb
    group: "1kG"
    shell:
        """
        cp {input} {params.out}
        """

rule convert_qced_snp_data_to_bfile_format:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/all/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/all/merged", ".bed", ".bim", ".fam")
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/all/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/all/merged"
    threads: 16
    resources:
        mem_mb = get_mem_mb
    group: "1kG"
    shell:
        """
        plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --make-bfile --silent --out {params.out_stem}
        """

rule create_pruned_ranges:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/merged", ".pgen", ".pvar.zst", ".psam")
    output:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out")
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged",
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 20
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --indep-pairwise {wildcards.window_size} 1 {params.r2} --out {params.out_stem}"

rule prune_snps:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/merged", ".pgen", ".pvar.zst", ".psam"),
        range_file = "results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged.prune.out"
    output:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged", ".pgen", ".psam", ".pvar.zst")
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/qc/{variant_set}/pruned/{window_size}_1_{r2}/merged"
    threads: 16
    resources:
        mem_mb = get_mem_mb
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --exclude {input.range_file} --make-pgen vzs --out {params.out_stem}"
