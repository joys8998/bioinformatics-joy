
rule process_intervals:
    input:
        ref="resources/genome.fa",
        intervals="resources/calling_regions.hg38.interval_list"
    params:
        bin_len=0,
        merging_rule="OVERLAPPING_ONLY",
        tmp_dir=config["tmp_dir"]
    output:
        interval_preprocessed="sandbox/calling_regions.preprocessed.interval_list"
    shell:
        """gatk PreprocessIntervals -L {input.intervals} \
        -R {input.ref} --bin-length {params.bin_len} \
        --interval-merging-rule {params.merging_rule} -O {output.interval_preprocessed} \
        --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*"""

rule collect_read_counts:
    input:
        bam="recal/{sample}-{unit}-{condition}.bam"
    params:
        intervals="sandbox/calling_regions.preprocessed.interval_list",
        merging_rule="OVERLAPPING_ONLY",
        tmp_dir=config["tmp_dir"]
    output:
        read_count="sandbox/{sample}-{unit}-{condition}.counts.hdf5"
    shell:
        """gatk CollectReadCounts -I {input.bam} \
        -L {params.intervals} \
        --interval-merging-rule {params.merging_rule} -O {output.read_count} \
        --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*"""

rule create_read_count_pon:
    input:
        counts=expand("sandbox/{u.sample}-{u.unit}-normal.counts.hdf5", u=units.itertuples())
    params:
        i=lambda wildcards, input: ['-I ' + d for d in input.counts],
        tmp_dir=config["tmp_dir"],
        mimp=5.0
    output:
        pon_hdf5="sandbox/cnv.pon.hdf5"
    shell:
        """gatk CreateReadCountPanelOfNormals --java-options '-Xmx6500m' {params.i} \
        --minimum-interval-median-percentile {params.mimp} \
        -O {output.pon_hdf5} --tmp-dir {params.tmp_dir}   && rm -r {params.tmp_dir}*      
        """

rule annotate_intervals:
    input:
        ref="resources/genome.fa"
    params:
        intervals="sandbox/calling_regions.preprocessed.interval_list",
        merging_rule="OVERLAPPING_ONLY",
        tmp_dir=config["tmp_dir"]
    output:
        annotated_intervals="sandbox/calling_regions_annotated_intervals.tsv"
    shell:
        """gatk AnnotateIntervals -R {input.ref} \
        -L {params.intervals} --interval-merging-rule {params.merging_rule} \
        -O {output.annotated_intervals} \
        --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*"""

rule denoise_read_count:
    input:
        annotated_intervals="sandbox/calling_regions_annotated_intervals.tsv",
        pon_hdf5="sandbox/cnv.pon.hdf5",
        read_count="sandbox/{sample}-{unit}-{condition}.counts.hdf5"
    params:
        tmp_dir=config["tmp_dir"]
    output:
        std_copy_ratios="sandbox/{sample}-{unit}-{condition}.standardizedCR.tsv",
        denoised_copy_ratios="sandbox/{sample}-{unit}-{condition}.denoisedCR.tsv"
    shell:
        """gatk DenoisereadCounts -I {input.read_count} \
        --denoised-copy-ratios {output.denoised_copy_ratios} --standardized-copy-ratios {output.std_copy_ratios} \
        --annotated-intervals {input.annotated_intervals} --tmp-dir {params.tmp_dir} --java-options '-Xmx12g' \
        && rm -r {params.tmp_dir}*"""

rule install_r_dependencies:
    output:
        "end_r.txt"
    conda:
        "../envs/r_plot.yml"
    shell:
        """Rscript scripts/install_R_packages.R && 'r installation ended' > {output}"""



rule plot_denoised_copy_ratios:
    input:
        end_r="end_r.txt",
        std_copy_ratios="sandbox/{sample}-{unit}-{condition}.standardizedCR.tsv",
        denoised_copy_ratios="sandbox/{sample}-{unit}-{condition}.denoisedCR.tsv",
        dict="resources/genome.dict"
    params:
        min_contig_len=46709983,
        tmp_dir=config["tmp_dir"],
        out_prefix="{sample}-{unit}-{condition}"
    output:
        dir=directory("sandbox/plots"),
        png="sandbox/plots/{sample}-{unit}-{condition}.denoised.png"
    shell:
        """rm {input.end_r} && gatk PlotDenoisedCopyRatios --standardized-copy-ratios {input.std_copy_ratios} \
        --denoised-copy-ratios {input.denoised_copy_ratios} \
        --sequence-dictionary {input.dict} --minimum-contig-length {params.min_contig_len} \
        --output {output.dir} --output-prefix {params.out_prefix} \
         --tmp-dir {params.tmp_dir} && rm -r {params.tmp_dir}*"""



