rule download_hg38_reference_sequence:
    output:
        ensure("resources/genome_reference/hg38.fa.zst", sha256 = "abf5670bf3de11a1e39176bc6596fb0597071e2be33ecae269855e502e1cdc53")
    params:
        url = "https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1"
    resources:
        runtime = 30
    localrule: True
    shell:
        """
        wget -O {output} {params.url}
        """

rule download_1kG_hg38_manifest:
    output:
        temp("resources/1kG/hg38/manifest.tsv")
    params:
        url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv"
    localrule: True
    shell:
        "wget -O {output} {params.url}"

rule process_1kG_hg38_manifest:
    input:
        "resources/1kG/hg38/manifest.tsv"
    output:
        "results/1kG/hg38/processed_manifest.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t', names = ['File', 'Byte', 'Checksum'])

        daf = daf[daf.File.str.match(r'.+\.vcf\.gz$')]

        daf = daf.assign(Chr=daf.File.str.extract(r'.+chr(\w+)\.filtered.+'))

        daf.to_csv(output[0], sep = '\t', index = False)

rule download_1kG_hg38_genotype_data:
    output:
        protected(ensure("resources/1kG/{assembly,hg38}/{chr}.vcf.gz", sha256 = get_vcf_sha256))
    params:
        chrX_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz",
        autosome_url = lambda w: f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{w.chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
    resources:
        runtime = 60
    retries: 3
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O resources/1kG/hg38/chrX.vcf.gz {params.chrX_url}")
        else:
            shell("wget -O resources/1kG/hg38/{wildcards.chr}.vcf.gz {params.autosome_url}")

rule download_1kG_hg38_sample_metadata:
     output:
        "resources/1kG/hg38/ped.txt"
     params:
         url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
     localrule: True
     shell:
         "wget -O resources/1kG/hg38/ped.txt {params.url}"

rule download_hg38_recombination_map:
    output:
        ensure("resources/1kG/hg38/genetic_map_hg38_withX.txt.gz", sha256 = "01ebe9f1e40a9b7e0bc96ae26543ca54702fc78b05c0d8726d8eab608c2c76c9")
    params:
        url = "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
    localrule: True
    shell:
        """
        wget -O {output} {params.url}
        """
