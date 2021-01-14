from google.cloud import storage

bucket_name = snakemake.params.bucket_name
source_blob_name = snakemake.params.ref_url
destination_file_name = snakemake.output.ref


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


if __name__ == "__main__":
    download_blob(bucket_name, source_blob_name, destination_file_name)