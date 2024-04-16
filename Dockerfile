# BSD 3-Clause License
#
# Copyright 2024 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Example command:
# $ git clone https://github.com/zhengzhenxian/ru.git
# $ cd ru
# $ docker build -f ./Dockerfile -t hkubal/ru:latest .
# $ docker run -it hkubal/ru:latest /opt/bin/ru --help

FROM continuumio/miniconda3:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        g++ \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install anaconda
RUN  conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda create -n ru python=3.9.0 -y

ENV PATH /opt/conda/envs/ru/bin:$PATH
ENV CONDA_DEFAULT_ENV ru

RUN /bin/bash -c "source activate ru" && \
    conda install -c conda-forge pypy3.6 -y && \
    conda install -c anaconda pigz -y && \
    conda install -c anaconda cffi=1.14.4 -y && \
    conda install -c conda-forge parallel=20191122 zstd -y && \
    conda install -c conda-forge -c bioconda samtools=1.15.1 -y && \
    conda install -c conda-forge -c bioconda whatshap=1.7 -y && \
    conda install -c conda-forge xz zlib bzip2 -y && \
    conda install -c conda-forge automake curl -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip && \
    echo "source activate ru" > ~/.bashrc

COPY . .

RUN cd /opt/bin && \
    make PREFIX=/opt/conda/envs/ru PYTHON=/opt/conda/envs/ru/bin/python && \
    rm -rf /opt/bin/samtools-*
