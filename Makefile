OS := $(shell uname)
ARCH := $(shell arch)

PYTHON ?= python3

all : libhts.a

SAMVER	=	1.15.1
LPVER	=	1.6
GCC	?=	gcc
GXX	?=	g++
PREFIX	?=	${CONDA_PREFIX}
LDFLAGS	=	-L ${PREFIX}/lib
CFLAGS	= -fpic -std=c99 -O3 -I ${PREFIX}/include -L ${PREFIX}/lib
CPPFLAGS	=	-std=c++11 -Wall -O3 -I ${PREFIX}/include -L ${PREFIX}/lib -Wl,-rpath=${PREFIX}/lib


samtools-$(SAMVER)/Makefile:
		curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
		tar -xjf samtools-${SAMVER}.tar.bz2; \
		rm samtools-${SAMVER}.tar.bz2

libhts.a: samtools-$(SAMVER)/Makefile
	# this is required only to add in -fpic so we can build python module
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd samtools-${SAMVER}/htslib-${SAMVER}; CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}" ./configure; make CFLAGS="${CFLAGS}" LDFLAGS="${LDFLAGS}"
	cp samtools-${SAMVER}/htslib-${SAMVER}/$@ $@