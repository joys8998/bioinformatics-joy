if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/{contig}.regions.bed"
        conda:
            "../envs/bedops.yml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"



rule mutect2:
    input:
        unpack(get_sample_bams),
        ref="resources/genome.fa",
        germ_res="resources/cont_res.vcf.gz",
        interval="interval-files/{interval}-scattered.interval_list"
    output:
        vcf="vcfs/{sample}.{unit}.{interval}.unfiltered.vcf",
        idx="vcfs/{sample}.{unit}.{interval}.unfiltered.vcf.idx",
        stats="vcfs/{sample}.{unit}.{interval}.unfiltered.vcf.stats",
        f1r2tar="vcfs/{sample}.{unit}.{interval}.f1r2.tar.gz"
    params:
        normal_name= ' ' if tumor_only else '-normal {sample}-normal',
        normal_input=lambda wildcards, input: '' if tumor_only else "-I " + input.normal,
        tumor_name='-tumor {sample}-tumor',
        pon="-pon resources/pon.vcf.gz" if config["use_pon"] else "",
        extra="",
        tmp_dir = config["tmp_dir"],
        interval="- L called/{interval}.regions.bed" if config["processing"].get("restrict-regions") else [],
        probability_threshold=0.002
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.tumor} \
            {params.normal_input} {params.normal_name} \
            -O {output.vcf} {params.tumor_name} \
            --germline-resource {input.germ_res} \
            --f1r2-tar-gz {output.f1r2tar} \
            {params.interval} \
            {params.pon} {params.extra}  -L {input.interval}\
            --active-probability-threshold {params.probability_threshold} \
            --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx75G -XX:ParallelGCThreads=10' \
            --tmp-dir {params.tmp_dir}
        """



rule orientation_bias:
    input:
        get_orientationbias_input
    output:
        "vcfs/{sample}.{unit}.read_orientation_model.tar.gz"
    params:
        i=lambda wildcards, input: ['-I ' + d for d in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk LearnReadOrientationModel {params.i} -O {output}
        """

rule pileup_summaries:
    input:
        bam="recal/{sample}-{unit}-{condition}.bam",
        cont_res="resources/cont_res.vcf.gz",
    output:
        "qc/{sample}-{unit}-{condition}_pileupsummaries.table"
    params:
        tmp_dir=config["tmp_dir"]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk GetPileupSummaries -I {input.bam} -V {input.cont_res} \
            -L {input.cont_res} -O {output} \
             --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx75G' \
            --tmp-dir {params.tmp_dir} 
        """

rule calculate_contamination:
    input:
        unpack(get_contamination_input)
    output:
        "qc/{sample}_{unit}_contamination.table"
    params:
        matched=lambda wildcards, input:'' if tumor_only else '-matched ' + input.normal
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk CalculateContamination -I {input.tumor}  \
            {params.matched} \
            -O {output}
        """

rule merge_vcfs:
    input:
        get_mergevcfs_input
    output:
        vcf="vcfs/{sample}.unfiltered.vcf"
    params:
        i=lambda wildcards, input: ['-I ' + vcf for vcf in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeVcfs {params.i} -O {output.vcf}
        """


rule merge_stats:
    input:
        get_mergestats_input
    output:
        stats="vcfs/{sample}.unfiltered.vcf.stats"
    params:
        i=lambda wildcards, input: ['-stats ' + s for s in input]
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk MergeMutectStats {params.i} -O {output.stats} 
        """


rule filter_calls:
    input:
        vcf="vcfs/{sample}.unfiltered.vcf",
        ref="resources/genome.fa",
        contamination="qc/{sample}_{unit}_contamination.table",
        stats="vcfs/{sample}.unfiltered.vcf.stats",
        f1r2model="vcfs/{sample}.read_orientation_model.tar.gz"
    output:
        vcf="vcfs/{sample}-{unit}.vcf",
        idx="vcfs/{sample}-{unit}.vcf.idx",
        intermediate=temp("vcfs/{sample}-{unit}.unselected.vcf"),
        inter_stats="vcfs/{sample}-{unit}.unselected.vcf.filteringStats.tsv",
        inter_idx=temp("vcfs/{sample}-{unit}.unselected.vcf.idx")
    conda:
        "../envs/gatk.yml"
    shell:
        """
        gatk FilterMutectCalls -V {input.vcf} -R {input.ref} \
            --contamination-table {input.contamination} \
            --stats {input.stats} \
            -ob-priors {input.f1r2model} \
            --min-reads-per-strand 1 \
            -O {output.intermediate}
        gatk SelectVariants -V {output.intermediate} -R {input.ref} -O {output.vcf} \
            --exclude-filtered -OVI 
        """


rule funcotator_datasource_downloader:
    output:
        directory("../DATA_SOURCES_DIR")
    conda:
        "../envs/gatk.yml"
    shell:
        "gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download"

rule enable_gnomAD:
    output:
        directory("../DATA_SOURCES_DIR/funcotator_dataSources.v1.6.20190124s/gnomAD_exome"),
        directory("../DATA_SOURCES_DIR/funcotator_dataSources.v1.6.20190124s/gnomAD_genome")
    shell:
        "tar -zxf ../DATA_SOURCES_DIR/funcotator_dataSources.v1.6.20190124s/gnomAD_exome.tar.gz && \
        tar -zxf ../DATA_SOURCES_DIR/funcotator_dataSources.v1.6.20190124s/gnomAD_genome.tar.gz"

rule annotation_funcotator:
    input:
        ref="resources/genome.fa",
        vcf="vcfs/{sample}-{unit}.vcf",
        data_sources_path="../DATA_SOURCES_DIR/funcotator_dataSources.v1.6.20190124s",
    params:
        ref_version="hg38",
        output_format="VCF",
        tmp_dir=config["tmp_dir"]
    conda:
        "../envs/gatk.yml"
    output:
        final="final/{sample}-{unit}-variants.funcotated.vcf"
    shell:
        """
        gatk Funcotator --variant {input.vcf} --reference {input.ref} --ref-version {params.ref_version} \
        --data-sources-path {input.data_sources_path} --output {output.final} --output-file-format {params.output_format} \
        --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' --tmp-dir {params.tmp_dir}
        """

