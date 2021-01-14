SAMPLES = ['A', 'B', 'C', 'D']
UNIT = 1
CONDITIONS = ['tumor', 'normal']
include: "rules/common.smk"

rule all:
	input:
		ref="resources/genome.fa"
		#seq_dict="resources/genome.dict",
		#cont_idx="resources/cont_res.vcf.gz.tbi",
		#pon_idx="resources/pon.vcf.gz.tbi",
		#genome_sa="resources/genome.fa.sa",
		#fai="resources/genome.fa.fai",
		#calling_regions="resources/calling_regions.hg38.interval_list",
		#interval="interval-files/0000-scattered.interval_list",
		#high_confidence="resources/high_confidence.vcf.gz",
		#bam_2=expand("recal/{sample}-{unit}-{condition}.bam", sample=SAMPLES, unit=UNIT, condition=CONDITIONS),
		#bai_2=expand("recal/{sample}-{unit}-{condition}.bam.bai",sample=SAMPLES,unit=UNIT,condition=CONDITIONS),
		#multi_qc="qc/multiqc.html",
		#vcf_merged=expand("vcfs/{sample}.{unit}.read_orientation_model.tar.gz", sample=SAMPLES, unit=UNIT),
		#contamination=expand("qc/{sample}_{unit}_contamination.table", sample=SAMPLES, unit=UNIT),
		#pileup=expand("qc/{sample}-{unit}-{condition}_pileupsummaries.table", sample=SAMPLES, unit=UNIT, condition=CONDITIONS),
		#vcf_filtered=expand("vcfs/{sample}-{unit}.vcf", sample=SAMPLES, unit=UNIT),
		#final=expand("final/{sample}-{unit}-variants.funcotated.vcf", sample=SAMPLES, unit=UNIT)
		# html=expand("qc/fastqc/{sample}.html", sample=SAMPLES),



# Modules
include: "rules/reference.smk"
include: "rules/preprocessing_mapping.smk"
include: "rules/qc.smk"
include: "rules/pon.smk"
include: "rules/calling.smk"


# # REMOVE SEC AND SUPP ALIGNMENTS
# # TRIM ADAPTER
# # MARK A/T AND C/G
# # REMOVE DUPLICATES
# rule path_seq_build_kmers:
#     input:
#         "../data-sn-tutorial/genome.fa"
#     output:
#         "../references/host_reference.hss.bfi"
#     shell:
#         "gatk PathSeqBuildKmers --reference {input} --output {output} --bloom-false-positive-probability 0.001 --kmer-mask 16 --kmer-size 31"
#
