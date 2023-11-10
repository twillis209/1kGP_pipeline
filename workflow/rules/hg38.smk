rule download_hg38_reference_sequence:
    output:
        ensure("resources/genome_reference/hg38.fa.zst", sha256 = "abf5670bf3de11a1e39176bc6596fb0597071e2be33ecae269855e502e1cdc53")
    resources:
        runtime = 30
    localrule: True
    shell:
        """
        wget -O {output} https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1
        """

rule download_1kG_hg38_manifest:
    output:
        temp("resources/1kG/hg38/manifest.tsv")
    localrule: True
    shell:
        "wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv"

rule process_1kG_hg38_manifest:
    input:
        "resources/1kG/hg38/manifest.tsv"
    output:
        "results/1kG/hg38/processed_manifest.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t', names = ['File', 'Byte', 'Checksum'])

        daf = daf[daf.File.str.match('.+\.vcf\.gz$')]

        daf = daf.assign(Chr=daf.File.str.extract('.+chr(\w+)\.filtered.+'))

        daf.to_csv(output[0], sep = '\t', index = False)

rule download_1kG_hg38_genotype_data:
    output:
        protected(ensure("resources/1kG/{assembly,hg38}/{chr}.vcf.gz", sha256 = get_vcf_sha256))
    resources:
        runtime = 60
    retries: 3
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O resources/1kG/hg38/chrX.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
        else:
            shell("wget -O resources/1kG/hg38/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

rule download_1kG_hg38_sample_metadata:
     output:
        "resources/1kG/hg38/ped.txt"
     localrule: True
     shell:
         "wget -O resources/1kG/hg38/ped.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule download_hg38_recombination_map:
    output:
        "resources/1kG/hg38/genetic_map_hg38_withX.txt.gz"
    localrule: True
    shell:
        """
        wget -O {output} https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
        """

rule write_out_per_chrom_hg38_recombination_map_files:
    input:
        "resources/1kG/hg38/genetic_map_hg38_withX.txt.gz"
    output:
        temp([f"resources/1kG/hg38/genetic_map_hg38/chr{x}.txt" for x in list(range(1,23))+['X']])
    params:
        root_dir = "resources/1kG/hg38/genetic_map_hg38"
    localrule: True
    threads: 1
    resources:
        mem_mb = get_mem_mb,
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
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged", ".bed", ".bim", ".fam"),
        map_files = [f"resources/1kG/hg38/genetic_map_hg38/chr{x}.txt" for x in list(range(1,23))+['X']]
    output:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged_with_cm", ".bed", ".bim", ".fam")
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged",
        out_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged_with_cm",
        map_pattern = "resources/1kG/hg38/genetic_map_hg38/chr@.txt"
    threads: 16
    resources:
        mem_mb = get_mem_mb,
        runtime = 5
    group: "1kG"
    conda: "../envs/1kGP_pipeline.yaml"
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --cm-map {params.map_pattern} --make-bed --out {params.out_stem}"
