rule download_hg19_variant_reference:
    output:
        "resources/genome_reference/hg19.tsv.gz"
    group: "1kG"
    shell:
        """
        wget -O {output} ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
        """

rule download_hg19_reference_sequence:
    output:
        ensure("resources/genome_reference/hg19.fa.zst", sha256 = "32f65df649ae46813bad00fee998542c7fd121aa9d01659e950ac307f2502693")
    params:
        compressed = "resources/genome_reference/hg19.fa.gz",
        uncompressed = "resources/genome_reference/hg19.fa"
    localrule: True
    shell:
        """
        wget -O {params.compressed} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
        gunzip {params.compressed}
        zstd --rm {params.uncompressed}
        """

rule download_1kG_hg19_genotype_data:
    output:
        protected(ensure("resources/1kG/{assembly,hg19}/{chr}.vcf.gz", sha256 = get_vcf_sha256))
    resources:
        runtime = 15
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz")
        elif wildcards.chr == 'chrY':
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz")
        else:
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

rule download_1kG_hg19_sample_metadata:
    output:
        panel = "resources/1kG/hg19/panel.txt",
        ped = "resources/1kG/hg19/ped.txt"
    localrule: True
    shell:
        """
        wget -O {output.panel} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        wget -O {output.ped} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped
        """
