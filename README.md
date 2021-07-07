# Representation Unification

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/clair3/README.html)

Contact: Ruibang Luo  

Email: rbluo@cs.hku.hk  

---

## Introduction

This document shows how Clair3 unifies the representation between the training materials and true variant set.

----

## Prerequisites

- Clair3 installed 
- GNU Parallel installed
- Pypy3 installed

----

## Input/Output

**Input:**

- BAM: Indexed BAM input
- REF: Indexed Reference input
- VCF: True variant VCF input
- BED: Confident BED regions input (optional)

**Ouput:**

- An VCF file of truth variants unified to the training materials