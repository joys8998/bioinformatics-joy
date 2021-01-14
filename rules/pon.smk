UNIT= 1
SAMPLES = ['A', 'B', 'C']

if config["use_pon"] and config["build_pon"]:
    rule mutect2_pon:
        input:
            bam="recal/{sample}-{unit}-normal.bam",
            ref="resources/genome.fa"
        params:
            tmp_dir=config["tmp_dir"]
        output:
            vcf="pon/{sample}-{unit}.pon.vcf.gz",
            idx="pon/{sample}-{unit}.pon.vcf.gz.tbi",
            stats="pon/{sample}-{unit}.pon.vcf.gz.stats"
        conda:
            "../envs/gatk.yml"
        shell:
            """
            gatk Mutect2 -I {input.bam} -R {input.ref} -O {output.vcf} \
            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
            -max-mnp-distance 0 --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*
            """

if config["use_pon"] and config["build_pon"]:
    rule gather_variants:
        input:
            vcfs=expand("pon/{sample}-{unit}.pon.vcf.gz", sample=SAMPLES, unit=UNIT),
            ref="resources/genome.fa",
            intervals="resources/calling_regions.hg38.interval_list"
        output:
            directory("pon/pon_db")
        params:
            vcfs=lambda wildcards, input: " -V ".join(input.vcfs),
            tmp_dir=config["tmp_dir"],
            batch_size=50
        shell:
            """
            gatk GenomicsDBImport -R {input.ref} --genomicsdb-workspace-path {output} \
                -V {params.vcfs} -L {input.intervals} --reader-threads 2 \
                --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx4g -Xms4g' \
                --batch-size {params.batch_size} --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*
            """

if config["use_pon"] and config["build_pon"]:
    rule create_pon:
        input:
            var="pon/pon_db",
            ref="resources/genome.fa",
            ger_res="resources/germ_res.vcf.gz",
        params:
            tmp_dir=config["tmp_dir"]
        output:
            vcf="resources/pon.vcf.gz",
            idx="resources/pon.vcf.gz.tbi"
        shell:
            """
            gatk CreateSomaticPanelOfNormals -R {input.ref} -V gendb://{input.var} -O {output.vcf} \
            --germline-resource {input.ger_res} --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*
            """