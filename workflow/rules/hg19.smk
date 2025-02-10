rule download_hg19_variant_reference:
    output:
        "resources/genome_reference/hg19.tsv.gz"
    params:
        url = "ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
    group: "1kG"
    shell:
        """
        wget -O {output} {params.url}
        """

rule download_hg19_reference_sequence:
    output:
        ensure("resources/genome_reference/hg19.fa.zst", sha256 = "32f65df649ae46813bad00fee998542c7fd121aa9d01659e950ac307f2502693")
    params:
        compressed = "resources/genome_reference/hg19.fa.gz",
        uncompressed = "resources/genome_reference/hg19.fa",
        url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    localrule: True
    shell:
        """
        wget -O {params.compressed} {params.url}
        gunzip {params.compressed}
        zstd --rm {params.uncompressed}
        """

rule download_1kG_hg19_genotype_data:
    output:
       protected(ensure("resources/1kG/{assembly,hg19}/{chr}.vcf.gz", sha256 = get_vcf_sha256))
    params:
        chrX_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz",
        chrY_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz",
        autosome_url = lambda w: f"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{w.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    resources:
        runtime = 15
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O {output} {params.chrX_url}")
        elif wildcards.chr == 'chrY':
            shell("wget -O {output} {params.chrY_url}")
        else:
            shell("wget -O {output} {params.autosome_url}")

rule download_1kG_hg19_sample_metadata:
    output:
        panel = "resources/1kG/hg19/panel.txt",
        ped = "resources/1kG/hg19/ped.txt"
    params:
        panel_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
        ped_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped"
    localrule: True
    shell:
        """
        wget -O {output.panel} {params.panel_url}
        wget -O {output.ped} {params.ped_url}
        """

rule download_hg19_recombination_map:
    output:
        ensure("resources/1kG/hg19/genetic_map_hg19_withX.txt.gz", sha256 = "7398f142f02815bdf563e125a0c34770f4024eb250ba57cb6ef81570c109fbbf")
    params:
        url = "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz"
    localrule: True
    shell:
        """
        wget -O {output} {params.url}
        """
