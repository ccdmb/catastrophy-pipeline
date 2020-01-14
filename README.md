# catastrophy-pipeline

Run the [catastrophy fungal trophy classifier](https://github.com/ccdmb/catastrophy) for many genomes.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)


## Introduction

This pipeline automates the process of running HMMER3 and CATAStrophy on many genomes.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
It comes with docker containers making installation trivial and results highly reproducible.

This documentation is a bit sparse right now.
We're all time poor.
If you are having trouble, please don't hesitate to raise an issue on github or email me.


## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run ccdmb/catastrophy-pipeline -profile test,<docker/singularity/conda>
```

iv. Start running your own analysis!

```bash
nextflow run ccdmb/catastrophy-pipeline -profile <docker/singularity/conda> --proteomes 'proteomes/*.fasta' --dbcan_version 8
```

## Parameters

| Parameter | default | description |
| :---      | :---    | :---        |
| `--help` | flag | Show help text and exit. |
| `--proteomes` | File or glob of files. | The proteins in fasta format that catastrophy should classify. Each file is treated as a separate organism. |
| `--dbcan_version` | 4,5,6,7 or 8 | The version of dbCAN to count CAZymes from. By default, will attempt to download the correct dbCAN database for this version and use that. |
| `--dbcan` | File | A copy of the dbCAN HMMER3 formatted database of CAZymes. If you use this option, you must also provide the `--dbcan_version` number. |
| `--dbcan_url` | URL | The url to download the dbCAN HMMER3 formatted database from. If you use this option, you must also provide the `--dbcan_version` number. |
| `--outdir` | `results` | The directory to write results to. |


I highly recommend that you use the `--dbcan_version` flag without `--dbcan` or `--dbcan_url`.
This will ensure that the correct model and database versions are used.

Nextflow can run tasks in parallel, I will add some documentation to this later.
Generally, it should do some of this automatically but more advanced things like using supercomputers or distributed cloud setups are more complex.
