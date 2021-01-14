import pandas as pd
from snakemake.utils import validate
import os
from snakemake.utils import min_version

min_version("5.9.1")

##### CONFIG FILE ####
###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit", "condition"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
validate(units, schema="../schemas/units.schema.yaml")


####### SET USEFUL VARIABLES ######
num_workers = config["num_workers"]

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    condition="|".join(units["condition"])


######## HELPER FUNCTIONS ########

# contigs in reference genome
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit, wildcards.condition), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit, condition):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit, condition), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{run}\tSM:{sample}-{condition}\tPL:{platform}'".format(
        sample=wildcards.sample,
        condition=wildcards.condition,
        run=units.loc[(wildcards.sample, wildcards.unit, wildcards.condition), "fq1"].split("/")[-1].split(".")[0],
        platform=units.loc[(wildcards.sample, wildcards.unit, wildcards.condition), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}-{condition}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}-{condition}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    files = dict()
    files['tumor'] = f"recal/{wildcards.sample}-{units.loc[wildcards.sample].unit[0]}-tumor.bam"
    if not config["tumor_only"]:
        files['normal'] = f"recal/{wildcards.sample}-{units.loc[wildcards.sample].unit[0]}-normal.bam"
    return files



def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default



def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}-{condition}.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}-{unit}-{condition}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def choose_alignment(wildcards):
    platform = units.loc[(wildcards.sample, wildcards.unit, wildcards.condition), "platform"]
    if platform == "PACBIO":
        return expand("mapped/{sample}-{unit}-{condition}.b2.bam", **wildcards)
    return expand("mapped/{sample}-{unit}-{condition}.bam", **wildcards)


def get_intervals():
    ints = []
    for i in range(num_workers):
        num_zeros = 4 - len(str(i))
        interval = '0' * num_zeros + str(i)
        ints.append(interval)
    return ints

def get_interval_files():
    ints = get_intervals()
    files = [i + '-scattered.interval_list' for i in ints]
    files = [os.path.join("interval-files", f) for f in files]
    return files

interval_files = get_interval_files()

def get_orientationbias_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.sample}.{wildcards.unit}.{i}.f1r2.tar.gz" for i in intervals]
    return files



def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket."""
    storage_client = storage.Client("rosy-petal-301413")
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)
    blob.download_to_filename(destination_file_name)

    print(
        "Blob {} downloaded to {}.".format(
            source_blob_name, destination_file_name
        )
    )

def isWGS(wildcards):
    seqtype = units.loc[(wildcards.sample, wildcards.unit, wildcards.condition), "seqtype"]
    return seqtype == "WGS"

tumor_only = config["tumor_only"]

def get_contamination_input(wildcards):
    out = {}
    out['tumor'] = f'qc/{wildcards.sample}-{units.loc[wildcards.sample].unit[0]}-tumor_pileupsummaries.table'
    if not tumor_only:
        out['normal'] = f'qc/{wildcards.sample}-{units.loc[wildcards.sample].unit[0]}-normal_pileupsummaries.table'
    return out

def get_mergevcfs_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.sample}.{i}.unfiltered.vcf" for i in intervals]
    return files


def there_WGS():
    seq_type = list(units["seq_type"])
    if "WGS" not in seq_type:
        return False

def get_mergestats_input(wildcards):
    intervals = get_intervals()
    files = [f"vcfs/{wildcards.sample}.{i}.unfiltered.vcf.stats" for i in intervals]
    return files

def get_tumor_bams(wildcards):
    return expand("recal/{sample}-{unit}-tumor.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit)

def get_normal_bam(wildcards):
    return expand("recal/{sample}-{unit}-normal.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit)