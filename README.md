# Representation Unification

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo  

Email: rbluo@cs.hku.hk  

---

## Introduction

This document shows how Clair3 unifies the representation between the training materials and true variant set.

----

## Contents
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
- [Usage](#usage)
- [Disclaimer](#disclaimer)

----

## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/ru). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/ru:latest \
  /opt/bin/ru \
  --bam_fn ${INPUT_DIR}/sample.bam \       ## use your bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \           ## use your reference file name here
  --truth_vcf_fn ${INPUT_DIR}/truth.vcf \  ## use your truth VCF file name here
  --threads ${THREADS} \                   ## maximum threads to be used
  --platform ${PLATFORM} \                 ## options: {ont, hifi, ilmn}
  --output_dir ${OUTPUT_DIR}               ## output path prefix 
```
### Option 2. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
# clone the repo
git clone https://github.com/zhengzhenxian/ru.git
cd ru

# build a docker image named hkubal/ru:latest
# might require docker authentication to build docker image 
docker build -f ./Dockerfile -t hkubal/ru:latest .

# run the docker image like option 1
docker run -it hkubal/ru:latest /opt/bin/ru --help
```


Check [Usage](#Usage) for more options.

----

## Usage

### General Usage

```bash
./ru \
  --bam_fn ${INPUT_DIR}/sample.bam \       ## use your bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \           ## use your reference file name here
  --truth_vcf_fn ${INPUT_DIR}/truth.vcf \  ## use your truth VCF file name here
  --threads ${THREADS} \                   ## maximum threads to be used
  --platform ${PLATFORM} \                 ## options: {ont, hifi, ilmn}
  --output_dir ${OUTPUT_DIR}               ## output path prefix 

## Final output file: ${OUTPUT_DIR}/unified.vcf.gz
```

### Options

**Required parameters:**

```bash

  -b BAM_FN, --bam_fn BAM_FN
                        BAM file input. The input file must be samtools indexed.
  -r REF_FN, --ref_fn REF_FN
                        FASTA reference file input. The input file must be samtools indexed.
  --truth_vcf_fn TRUTH_VCF_FN
                        Truth VCF file input.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory.
  -t THREADS, --threads THREADS
                        Max threads to be used.
  -p PLATFORM, --platform PLATFORM
                        Select the sequencing platform of the input. Possible options: {ont, hifi,
                        ilmn}.

```

**Miscellaneous parameters:**

```bash

  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs.
                        Default: all contigs will be processed.
  --bed_fn BED_FN       Path to a BED file. Execute RU only in the provided BED regions.
  --region REGION       A region to be processed. Format: `ctg_name:start-end` (start is 1-based).
  --min_af MIN_AF       Minimal AF required for a variant to be called. Default: 0.08.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename. Default: output.
  --remove_intermediate_dir
                        Remove intermediate directory before finishing to save disk space.
  --include_all_ctgs    Execute RU on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}.
  -d, --dry_run         Print the commands that will be ran.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --whatshap WHATSHAP   Absolute path of whatshap, whatshap >= 1.0 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --disable_phasing     Disable phasing with whatshap.

```

----

## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
