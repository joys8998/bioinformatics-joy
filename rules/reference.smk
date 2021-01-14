# GET GENOME FROM ENSEMBL
rule get_genome:
    output:
        ref="resources/genome.fa"
    log:
        "logs/get-genome.log"
    params:
        ref_url= "references/hg38/v0/GRCh38.primary_assembly.genome.fa" if config["ref"]["build"] == "GRCh38" else "hg19/v0/Homo_sapiens_assembly19.fasta",
        bucket_name="genomics-public-data" if config["ref"]["build"] == "GRCh38" else "gcp-public-data--broad-references"
    cache: True
    conda:
        "../envs/gcs.yml"
    script:
        "../scripts/download_from_storage.py"
        #download_blob(params.bucket_name, params.ref_url, output.ref)


# CREATE FAI FILE
checkpoint genome_faidx:
    input:
        "resources/genome.fa"
    output:
        "resources/genome.fa.fai"
    log:
        "logs/genome-faidx.log"
    cache: True
    wrapper:
        "0.59.2/bio/samtools/faidx"

# CREATE REF DICTIONARY
rule genome_dict:
    input:
        "resources/genome.fa"
    output:
        "resources/genome.dict"
    log:
        "logs/picard/create_dict.log"
    cache: True
    wrapper:
        "0.67.0/bio/picard/createsequencedictionary"



# GET ALL KNOWN VARIANTS
if config["ref"]["build"] == "GRCh38":
    rule get_known_variation_hg38:
        output:
            "resources/high_confidence.vcf.gz",
            "resources/high_confidence.vcf.gz.tbi",
            "resources/dbsnp.vcf",
            "resources/dbsnp.vcf.idx",
            "resources/known_indels.vcf.gz",
            "resources/known_indels.vcf.gz.tbi",
        log:
            "logs/get-known-variants.log"
        params:
            "references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
            "references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
            "references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
            "references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
            "references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz",
            "references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
        cache: True
        run:
            for i, ele in enumerate(params):
                download_blob("genomics-public-data",params[i],output[i])
else:
    rule get_known_variation_hg19:
        output:
            "resources/high_confidence.vcf.gz",
            "resources/high_confidence.vcf.gz.tbi",
            "resources/dbsnp.vcf",
            "resources/dbsnp.vcf.idx",
        log:
            "logs/get-known-variants.log"
        params:
            "hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz",
            "hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz.tbi",
            "hg19/v0/Homo_sapiens_assembly19.dbsnp138.vcf",
            "hg19/v0/Homo_sapiens_assembly19.dbsnp138.vcf.idx",
        cache: True
        run:
            for i, ele in enumerate(params):
                download_blob("gcp-public-data--broad-references",params[i],output[i])

# GET WES INTERVALS
rule get_genome_intervals:
    output:
        calling_regions="resources/calling_regions.hg38.interval_list"
    params:
        link="resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list" if there_WGS else "hg38/v0/HybSelOligos/whole_exome_illumina_coding_v1/whole_exome_illumina_coding_v1.Homo_sapiens_assembly38.targets.interval_list",
        bucket_name="genomics-public-data" if there_WGS else "gcp-public-data--broad-references"
    run:
        download_blob(params.bucket_name, params.link, output.calling_regions)

# SPLIT THE INTERVALS BY THE NUMBERS OF WORKERS OF MUTECT2
rule split_intervals:
    input:
        ref="resources/genome.fa",
        intervals="resources/calling_regions.hg38.interval_list"
    output:
        interval_files
    params:
        N=num_workers,
        d="interval-files"
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk SplitIntervals -R {input.ref} -L {input.intervals} \
            --scatter-count {params.N} -O {params.d} \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
            --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
        """


# BWA indexing
rule bwa_index:
    input:
        "resources/genome.fa"
    output:
        multiext("resources/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    resources:
        mem_mb=369000
    cache: True
    wrapper:
        "0.59.2/bio/bwa/index"



rule download_germline_and_contamination_res:
    output:
        "resources/germ_res.vcf.gz",
        "resources/germ_res.vcf.gz.tbi",
        "resources/cont_res.vcf.gz",
        "resources/cont_res.vcf.gz.tbi"
    params:
        "somatic-hg38/af-only-gnomad.hg38.vcf.gz" if config["ref"]["build"] == "GRCh38" else "somatic-b37/af-only-gnomad.raw.sites.vcf",
        "somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" if config["ref"]["build"] == "GRCh38" else "somatic-b37/af-only-gnomad.raw.sites.vcf.tbi",
        "somatic-hg38/small_exac_common_3.hg38.vcf.gz" if config["ref"]["build"] == "GRCh38" else "somatic-b37/small_exac_common_3.vcf",
        "somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi" if config["ref"]["build"] == "GRCh38" else "somatic-b37/small_exac_common_3.vcf.tbi",
    run:
        for i, ele in enumerate(params):
            download_blob("gatk-best-practices",params[i],output[i])


if config["use_pon"] and config["build_pon"] == False:
    rule download_pon:
        output:
            "resources/pon.vcf.gz",
            "resources/pon.vcf.gz.tbi"
        params:
            "somatic-hg38/1000g_pon.hg38.vcf.gz" if config["ref"]["build"] == "GRCh38" else "somatic-b37/Mutect2-exome-panel.vcf",
            "somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" if config["ref"]["build"] == "GRCh38" else "somatic-b37/Mutect2-exome-panel.vcf.idx",
        run:
            for i, ele in enumerate(params):
                download_blob("gatk-best-practices",params[i],output[i])
