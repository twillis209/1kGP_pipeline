import polars as pl

def get_variant_set_filter_flags(wildcards, input):
    """
    Returns a string of plink2 flags for variant set filtering
    """
    variant_set_options = wildcards.variant_set.split("_and_")

    plink_flags = ""

    if 'all' in variant_set_options and len(variant_set_options) == 1:
        return plink_flags
    elif 'all' in variant_set_options:
        raise ValueError("Invalid variant set options, cannot specify 'all' with other options")

    # if 'sans_long_range_ld' or 'sans_at_gc' in variant_set_options:
    #     plink_flags += f" --exclude {input.snps_to_exclude}"

    if 'sans_long_range_ld' in variant_set_options:
        plink_flags += f" --exclude bed1 {input.long_range_ld}"

    if 'sans_mhc' in variant_set_options:
        plink_flags += f" --exclude bed1 {input.mhc}"

    if 'sans_pars' in variant_set_options:
        plink_flags += f" --not-chr PAR1 PAR2"

    return plink_flags

rule make_1kG_sex_file:
    input:
        "resources/1kG/{assembly}/ped.txt",
    output:
        "results/1kG/{assembly}/sex.tsv"
    localrule: True
    run:
        ped = pd.read_csv(input[0], sep = ' ', header = 0)

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
    log:
        "results/1kG/{assembly}/{chr}.log"
    params:
        out = subpath(output[0], strip_suffix = '.pgen'),
        id_format = r"@:#:\$r:\$a",
        max_allele_len = 20
    threads: 8
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input.vcf} --make-pgen vzs --fa {input.ref} --ref-from-fa 'force' --out {params.out} --set-all-var-ids {params.id_format} --max-alleles 2 --new-id-max-allele-len {params.max_allele_len} truncate --update-sex {input.sex} --split-par '{wildcards.assembly}'"

rule make_1kG_sample_files:
     input:
         ped = "resources/1kG/{assembly}/ped.txt",
         # No hg38 panel file required
         panel = "resources/1kG/hg19/panel.txt",
     params:
         out_dir = subpath(output[0], parent = True)
     output:
         expand("results/1kG/{{assembly}}/{{relatedness}}/{ancestry}.samples", ancestry = ["eur", "afr", "amr", "eas", "sas", "all"])
     localrule: True
     conda: env_path("r.yaml")
     script: "../scripts/write_1kG_sample_files.R"

rule get_ancestry_specific_samples:
     input:
        rules.vcf_to_pgen.output,
        sample_file = "results/1kG/{assembly}/{relatedness}/{ancestry}.samples"
     output:
        temp(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{chr}", ".pgen", ".pvar.zst", ".psam"))
     log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{chr}.log"
     params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen')
     threads: 8
     group: "1kG"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --keep {input.sample_file} --make-pgen vzs --out {params.out_stem}"

rule filter_for_snps_and_maf:
    input:
        rules.get_ancestry_specific_samples.output
    output:
        temp(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/{chr}", ".pgen", ".pvar.zst", ".psam"))
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/{chr}.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen'),
        snp_flag = lambda w: "--snps-only 'just-acgt'" if w.variant_type == "snps_only" else "",
        maf_flag = lambda w: f"--maf 0.{w.maf}" if w.maf != "0" else ""
    threads: 8
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs {params.snp_flag} --rm-dup 'force-first' {params.maf_flag} --make-pgen vzs --out {params.out_stem}"

rule merge_pgen_files:
    input:
        expand("results/1kG/{{assembly}}/{{relatedness}}/{{ancestry}}/{{variant_type}}/{{maf}}/{chr}.{ext}", chr = [f"chr{x}" for x in range(1,23)]+["chrX"], ext = ["pgen", "pvar.zst", "psam"])
    output:
        pfiles = protected(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/merged", ".pgen", ".pvar.zst", ".psam")),
        pmerge_file = "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/pmerge.txt"
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/merged.log"
    params:
        in_dir = subpath(input[0], parent = True),
        out_stem = subpath(output[0], strip_suffix = '.pgen'),
        max_allele_ct = 2
    threads: 16
    group: "1kG"
    shell: """
        for i in {{1..22}}; do
        echo "chr$i" >>{output.pmerge_file}
        done

        echo "chrX" >>{output.pmerge_file}

        plink2 --memory {resources.mem_mb} --threads {threads} --pmerge-list {output.pmerge_file} pfile-vzs --pmerge-list-dir {params.in_dir} --merge-max-allele-ct {params.max_allele_ct} --pmerge-output-vzs --out {params.out_stem}
        """

# NB: supporting SNPs only atm
rule pgen_to_hap_and_legend:
    input:
        rules.vcf_to_pgen.output
    output:
        temp(multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type,snps_only}/{maf}/{chr}", ".haps", ".legend", ".sample"))
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/{chr}.log"
    params:
        stem = subpath(input[0], strip_suffix = '.pgen')
    threads: 16
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --export hapslegend --out {params.stem}"

rule concatenate_legend_files:
    input:
        expand("results/1kG/{{assembly}}/{{ancestry}}/{{variant_type}}/{{maf}}/{chr}.legend", chr = [f"chr{x}" for x in range(1,23)])
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/combined.legend.gz"
    params:
        uncompressed_output = subpath(output[0], strip_suffix = '.gz')
    localrule: True
    shell:
        """
        for x in {input}; do
        tail -n +2 $x >>{params.uncompressed_output}
        done

        gzip {params.uncompressed_output}
        """

rule qc:
     input:
         rules.merge_pgen_files.output.pfiles
     output:
         multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type,snps_only}/{maf}/qc/merged", ".pgen", ".pvar.zst", ".psam")
     log:
         "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/merged.log"
     params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen'),
        geno = 0.01,
        mind = 0.01,
        hwe = 1e-50
     threads: 16
     resources:
        runtime = 10
     group: "1kG"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --geno {params.geno} --mind {params.mind} --hwe {params.hwe} --make-pgen vzs --out {params.out_stem}"

rule decompress_pvar_for_variant_screens:
    input:
        rules.qc.output[1]
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type,snps_only}/{maf}/qc/merged.pvar"
    localrule: True
    shell:
        "zstdcat {input} | grep -v '^#' >{output}"

rule identify_at_gc_snps:
    input:
        rules.decompress_pvar_for_variant_screens.output
    output:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/at_gc_snps.txt"
    threads: 4
    run:
        daf = pl.read_csv(input[0], separator = '\t', columns = [2, 3, 4], has_header = False, new_columns = ['ID', 'REF', 'ALT'])

        ambiguous = daf.filter(
            ((pl.col("REF") == "A") & (pl.col("ALT") == "T")) |
            ((pl.col("REF") == "T") & (pl.col("ALT") == "A")) |
            ((pl.col("REF") == "C") & (pl.col("ALT") == "G")) |
            ((pl.col("REF") == "G") & (pl.col("ALT") == "C"))
        )

        ambiguous.select(pl.col("ID")).write_csv(output[0], separator = '\t')

# TODO not sure how to easily combine bed regions with individual SNPs
# rule identify_snps_to_exclude:
#     input:
#         at_gc_snps = rules.identify_at_gc_snps.output,
#         long_range_ld = "resources/1kG/{assembly}/long_range_ld_regions.bed"
#     output:
#         "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/snps_to_exclude.txt"
#     localrule: True
#     run:
#         shell("touch {output}")

#         if 'sans_mhc' in wildcards.variant_set:
#             shell("tail -n +2 {input.mhc_snps} >> {output}")
#         if 'sans_at_gc' in wildcards.variant_set:
#             shell("tail -n +2 {input.at_gc_snps} >> {output}")

rule filter_variant_set:
    input:
        pfiles = rules.qc.output,
        long_range_ld = "resources/1kG/{assembly}/long_range_ld_regions.bed",
        mhc = "resources/1kG/{assembly}/mhc.bed"
    output:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged", ".pgen", ".pvar.zst", ".psam")
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen'),
        filter_flags = lambda w, input: get_variant_set_filter_flags(w, input)
    threads: 16
    resources:
        runtime = 10
    group: "1kG"
    run:
        shell("plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --make-pgen vzs --out {params.out_stem} {params.filter_flags}")

rule write_out_qced_data_to_bed_format:
    input:
        rules.filter_variant_set.output
    output:
        multiext("results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged", ".bed", ".bim", ".fam")
    log:
        "results/1kG/{assembly}/{relatedness}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.bed')
    threads: 16
    group: "1kG"
    shell:
        """
        plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --make-bfile --silent --out {params.out_stem}
        """
