SAMPLES = ['A', 'B', 'C', 'D']
UNITS = 1
CONDITIONS = ['tumor', 'normal']

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}-{unit}-{condition}.html",
        zip="qc/fastqc/{sample}-{unit}-{condition}.zip"
    wrapper:
        "0.59.2/bio/fastqc"



rule fastqc_trimmed:
    input:
        get_trimmed_reads
    output:
        html="qc/trimmed/{sample}-{unit}-{condition}.html",
        zip="qc/trimmed/{sample}-{unit}-{condition}.zip"
    wrapper:
        "0.59.2/bio/fastqc"

rule samtools_stats:
    input:
        "recal/{sample}-{unit}-{condition}.bam"
    output:
        "qc/samtools-stats/{sample}-{unit}-{condition}.txt"
    log:
        "logs/samtools-stats/{sample}-{unit}-{condition}.log"
    wrapper:
        "0.59.2/bio/samtools/stats"

rule multiqc:
    input:
        expand(["qc/samtools-stats/{u.sample}-{u.unit}-{u.condition}.txt",
                "qc/fastqc/{u.sample}-{u.unit}-{u.condition}.zip",
                "qc/dedup/{u.sample}-{u.unit}-{u.condition}.metrics.txt",
                "qc/trimmed/{u.sample}-{u.unit}-{u.condition}.zip"],
            u=units.itertuples()),
    output:
        "qc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"


