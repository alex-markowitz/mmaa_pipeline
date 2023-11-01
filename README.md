# MMAA Pipeline

Execute MMAA Pipeline for DNA methylation-based classification of central nervous system tumours.

## Author
- **Alexander Markowitz**

## Publication Reference
- Capper et al., Nature, 2018. [PubMed 29539639](https://pubmed.ncbi.nlm.nih.gov/29539639/)

## Overview

This script processes and classifies DNA methylation data for central nervous system tumours. It reads in a sample sheet, performs the required pre-processing on the samples, classifies the samples, performs CNV analysis, MGMT classification, and finally, generates a quality control report for each sample.

## Dependencies

- `optparse` library for R
- External utility scripts and data files located in a specified utility directory.

## Docker Usage

This pipeline is also available as a Docker container. You can pull and run the container using the following commands:

```bash
docker pull alexmarkowitz/mmaa_pipeline
docker run -v /path/to/your/data:/data alexmarkowitz/mmaa_pipeline Rscript /path/in/container/<script_name>.R -s /data/<sample_sheet.csv> [other options]
```

Visit [Docker Hub](https://hub.docker.com/r/alexmarkowitz/mmaa_pipeline) for more details on the Docker image.

## Script Usage

```bash
Rscript <script_name>.R -s <sample_sheet.csv> [-u <util_dir>] [-C] [-R] [-o <output_dir>]
```

### Parameters

| Short Flag | Long Flag      | Description                                                                                                                                      |
|------------|----------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| -s         | --sample-sheet | REQUIRED: CSV file with the sample sheet. Must contain 5 columns: sentrix barcode, sample ID, sample type (FFPE or Frozen), IDAT Red and IDAT Grn. |
| -u         | --util-dir     | Path to utility directory containing R scripts and data files. Default: `/usr/local/util_files/`                                                  |
| -C         | --CNV          | Flag to run CNV evaluation. Default: `FALSE`                                                                                                     |
| -R         | --RData        | Flag to save .RData session image at the end of the analysis. Default: `FALSE`                                                                   |
| -o         | --output-dir   | Directory name for HTML output files. Default: `/docker_scratch`                                                                                 |
