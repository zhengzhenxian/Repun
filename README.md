# Repun: An accurate small variant representation unification method for multiple sequencing platforms


[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Zhenxian Zheng, Ruibang Luo  

Email: {zxzheng,rbluo}@cs.hku.hk  

---

## Introduction

Ensuring a unified variant representation aligning the sequencing data is critical for downstream analysis as variant representation may differ across platforms and sequencing conditions. Current approaches typically treat variant unification as a post-step following variant calling and are incapable of measuring the correct variant representation from the outset. Aligning variant representations with the alignment before variant calling has benefits like providing reliable training labels for deep learning-based variant caller model training and enabling direct assessment of alignment quality. However, it also poses challenges due to the large number of candidates to handle. Here, we present Repun, a haplotype-aware variant-alignment unification algorithm that harmonizes the variant representation between provided variants and alignments in different sequencing platforms. Repun leverages phasing to facilitate equivalent haplotype matches between variants and alignments. Our approach reduced the comparisons between variant haplotypes and candidate haplotypes by utilizing haplotypes with read evidence to speed up the unification process. Repun achieved >99.99% precision and >99.5% recall through extensive evaluations of various GIAB samples encompassing three sequencing platforms: ONT, PacBio, and Illumina.



----

## Contents
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Docker Dockerfile](#option-2-docker-dockerfile)
- [Usage](#usage)
- [Disclaimer](#disclaimer)

----

## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/repun). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker n -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/repun:latest \
  /opt/bin/repun \
  --bam_fn ${INPUT_DIR}/sample.bam \       ## use your bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \           ## use your reference file name here
  --truth_vcf_fn ${INPUT_DIR}/tth.vcf \  ## use your VCF file name here
  --threads ${THREADS} \                   ## maximum threads to be used
  --platform ${PLATFORM} \                 ## options: {ont, hifi, ilmn}
  --output_dir ${OUTPUT_DIR}               ## output path prefix 
```
### Option 2. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
# clone the repo
git clone https://github.com/zhengzhenxian/Repun.git
cd Repun

# build a docker image named hkubal/repun:latest
# might require docker authentication to build docker image 
docker build -f ./Dockerfile -t hkubal/repun:latest .

# run the docker image like option 1
docker run -it hkubal/repun:latest /opt/bin/repun --help
```


Check [Usage](#Usage) for more options.

----

## Usage

### General Usage

```bash
./repun \
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
  --bed_fn BED_FN       Path to a BED file. Execute Repun only in the provided BED regions.
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
  --include_all_ctgs    Execute Repun on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}.
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
