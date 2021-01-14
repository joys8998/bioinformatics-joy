# TRIMMOMATIC
rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}-{unit}-{condition}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    threads:
        32
    log:
        "logs/trimmomatic/{sample}-{unit}-{condition}.log"
    wrapper:
        "0.59.2/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}-{condition}.1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}-{condition}.2.fastq.gz"),
        r1_unpaired=temp("trimmed/{sample}-{unit}-{condition}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("trimmed/{sample}-{unit}-{condition}.2.unpaired.fastq.gz"),
        trimlog="trimmed/{sample}-{unit}-{condition}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    threads:
        32
    log:
        "logs/trimmomatic/{sample}-{unit}-{condition}.log"
    wrapper:
        "0.59.2/bio/trimmomatic/pe"


# ALIGN SHORT READS
rule bowtie2_build:
    input:
        reference="resources/genome.fa"
    output:
        multiext(
            "resources/genome",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    threads: 8
    wrapper:
        "0.68.0/bio/bowtie2/build"


rule bowtie2:
    input:
        sample=get_trimmed_reads,
        genome="resources/genome.1.bt2"
    output:
        temp("mapped/{sample}-{unit}-{condition}.b2.bam")
    log:
        "logs/bowtie2/{sample}-{unit}-{condition}.b2.log"
    params:
        index="resources/genome",
        extra=""  #
    threads: 8  #
    wrapper:
        "0.68.0/bio/bowtie2/align"




# ALIGN NORMAL LENGTH READS AND SORT THEM USING SAMTOOLS
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output
    output:
        temp("mapped/{sample}-{unit}-{condition}.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}-{condition}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.59.2/bio/bwa/mem"



rule mark_duplicates:
    input:
        choose_alignment
    output:
        bam=temp("dedup/{sample}-{unit}-{condition}.bam"),
        metrics="qc/dedup/{sample}-{unit}-{condition}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}-{condition}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.59.2/bio/picard/markduplicates"






rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fa",
        known_sites=["resources/high_confidence.vcf.gz",
            "resources/dbsnp.vcf",
            "resources/known_indels.vcf.gz"]
    output:
        bam=protected("recal/{sample}-{unit}-{condition}.bam"),
        bai="recal/{sample}-{unit}-{condition}.bai",
        md5="recal/{sample}-{unit}-{condition}.bam.md5",
        recal="qc/{sample}-{unit}-{condition}.recal_data.table"
    params:
        ks=lambda wildcards, input: ["--known-sites " + s for s in input.known_sites]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} -O {output.recal} \
            {params.ks} 
        gatk ApplyBQSR -I {input.bam} -R {input.ref} -O {output.bam} -bqsr {output.recal} \
            --add-output-sam-program-record \
            --create-output-bam-md5
        """

rule samtools_index:
    input:
        "recal/{sample}-{unit}-{condition}.bam"
    output:
        "recal/{sample}-{unit}-{condition}.bam.bai"
    wrapper:
        "0.59.2/bio/samtools/index"