#!/usr/bin/env python

# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
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

import os
import sys
import argparse
import shlex
import subprocess

from collections import defaultdict, namedtuple
from argparse import SUPPRESS
try:
    from packaging.version import parse as version_parse
except ModuleNotFoundError:
    from distutils.version import LooseVersion as version_parse

from time import time

import shared.param as param
from shared.interval_tree import bed_tree_from
from shared.utils import file_path_from, folder_path_from, subprocess_popen, str2bool, legal_range_from, log_error, log_warning

major_contigs = {"chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]}.union(
    {str(a) for a in list(range(1, 23)) + ["X", "Y"]})
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]

file_directory = os.path.dirname(os.path.realpath(__file__))
main_entry = os.path.join(file_directory, "ru.py")
MAX_STEP = 10

OutputPath = namedtuple('OutputPath', [
    'log_path',
    'tmp_file_path',
    'split_bed_path',
    'candidates_path',
    'var_path',
    'phased_vcf_path',
    'phased_bam_path',
    'vcf_output_path',
])


class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()

def logging(str):
    if args.tee is None:
        print(str)
    else:
        args.tee.stdin.write(bytes(str + '\n', encoding='utf8'))

def create_output_folder(args):
    # create temp file folder
    args.output_dir = folder_path_from(args.output_dir, create_not_found=True)
    log_path = folder_path_from(os.path.join(args.output_dir, 'logs'), create_not_found=True)
    tmp_file_path = folder_path_from(os.path.join(args.output_dir, 'tmp'), create_not_found=True)
    split_bed_path = folder_path_from(os.path.join(tmp_file_path, 'split_beds'), create_not_found=True)
    candidates_path = folder_path_from(os.path.join(tmp_file_path, 'candidates'), create_not_found=True)
    var_path = folder_path_from(os.path.join(tmp_file_path, 'var'), create_not_found=True)
    phased_vcf_path = folder_path_from(os.path.join(tmp_file_path, 'phased_vcf'), create_not_found=True)
    phased_bam_path = folder_path_from(os.path.join(tmp_file_path, 'phased_bam'), create_not_found=True)
    vcf_output_path = folder_path_from(os.path.join(tmp_file_path, 'vcf_output'), create_not_found=True)

    output_path = OutputPath(log_path=log_path,
                             tmp_file_path=tmp_file_path,
                             split_bed_path=split_bed_path,
                             candidates_path=candidates_path,
                             var_path=var_path,
                             phased_vcf_path=phased_vcf_path,
                             phased_bam_path=phased_bam_path,
                             vcf_output_path=vcf_output_path)
    return output_path

def check_version(tool, pos=None, is_pypy=False):
    try:
        if is_pypy:
            proc = subprocess.run("{} -c 'import sys; print (sys.version)'".format(tool), stdout=subprocess.PIPE,
                                  shell=True)
        else:
            proc = subprocess.run([tool, "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        first_line = proc.stdout.decode().split("\n", 1)[0]
        version = first_line.split()[pos]
        version = version_parse(version)
    except Exception:
        return None

    return version

def check_skip_steps_legal(args):
    skip_steps = args.skip_steps
    skip_steps_list = skip_steps.rstrip().split(",")
    if len(skip_steps_list) == 0:
        sys.exit(log_error("[ERROR] --skip_steps option provided but no skip steps index found"))
    for step in skip_steps_list:
        if int(step) < 1 or int(step) > MAX_STEP:
            sys.exit(log_error("[ERROR] --skip_steps option provided but contains invalid skip steps index, should be 1-index"))

def check_python_path():
    python_path = subprocess.run("which python", stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
    sys.exit(log_error("[ERROR] Current python execution path: {}".format(python_path)))

def check_python_version(python):
    python_path = subprocess.run("{} --version".format(python), stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
    return python_path.split(' ')[1]


def check_tools_version(args):

    required_tool_version = {
        'python': version_parse('3.9.0'),
        'pypy': version_parse('3.6'),
        'samtools': version_parse('1.10'),
        'whatshap': version_parse('1.0'),
        'parallel': version_parse('20191122'),
    }

    tool_version = {
        'python': version_parse(check_python_version(args.python)),
        'pypy': check_version(tool=args.pypy, pos=0, is_pypy=True),
        'samtools': check_version(tool=args.samtools, pos=1),
        'parallel': check_version(tool=args.parallel, pos=2),
    }

    for tool, version in tool_version.items():
        required_version = required_tool_version[tool]
        if version is None:
            logging(log_error("[ERROR] {} not found, please check if you are in the conda virtual environment".format(tool)))
            check_python_path()
        elif version < required_version:
            logging(log_error("[ERROR] Tool version not match, please check if you are in the conda virtual environment"))
            logging(' '.join([str(item).ljust(10) for item in ["Tool", "Version", "Required"]]))
            error_info = ' '.join([str(item).ljust(10) for item in [tool, version, '>=' + str(required_version)]])
            logging(error_info)
            check_python_path()
    return


def check_contig_in_bam(bam_fn, sorted_contig_list, samtools, allow_none=False):
    if allow_none and bam_fn is None:
        return sorted_contig_list, True
    bai_process = subprocess_popen(shlex.split("{} idxstats {}".format(samtools, bam_fn)))
    contig_with_read_support_set = set()
    for row_id, row in enumerate(bai_process.stdout):
        row = row.split('\t')
        if len(row) != 4:
            continue
        contig_name, contig_length, mapped_reads, unmapped_reads = row
        if contig_name not in sorted_contig_list:
            continue
        if int(mapped_reads) > 0:
            contig_with_read_support_set.add(contig_name)
    for contig_name in sorted_contig_list:
        if contig_name not in contig_with_read_support_set:
            logging(log_warning(
                "[WARNING] Contig name {} provided but no mapped reads found in BAM, skip!".format(contig_name)))
    filtered_sorted_contig_list = [item for item in sorted_contig_list if item in contig_with_read_support_set]

    found_contig = True
    if len(filtered_sorted_contig_list) == 0:
        found_contig = False
        logging(log_warning(
            "[WARNING] No mapped reads found in BAM for provided contigs set {}".format(' '.join(sorted_contig_list))))

    return filtered_sorted_contig_list, found_contig


def check_threads(args):
    threads = args.threads
    #sched_getaffinity is not exist in pypy
    try:
        sched_getaffinity_list = list(os.sched_getaffinity(0))
        num_cpus = len(sched_getaffinity_list)
    except:
        num_cpus = int(subprocess.run(args.python + " -c \"import os; print(len(os.sched_getaffinity(0)))\"", \
                                      stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip())

    if threads > num_cpus:
        logging(log_warning(
            '[WARNING] Threads setting {} is larger than the number of available threads {} in the system,'.format(
                threads, num_cpus)))
        logging(log_warning('Set --threads={} for better parallelism.'.format(num_cpus)))
        args.threads = num_cpus
    return args

def split_extend_vcf(truth_vcf_fn, output_fn):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (truth_vcf_fn)))

    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split(maxsplit=3)
        ctg_name = columns[0]

        center_pos = int(columns[1])
        ctg_start, ctg_end = center_pos - 1, center_pos
        if ctg_start < 0:
            sys.exit(
                log_error("[ERROR] Invalid VCF input at the {}-th row {} {}".format(row_id + 1, ctg_name, center_pos)))
        if ctg_start - expand_region_size < 0:
            continue
        expand_ctg_start = ctg_start - expand_region_size
        expand_ctg_end = ctg_end + expand_region_size

        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()

    know_vcf_contig_set = set(list(output_ctg_dict.keys()))

    return know_vcf_contig_set


def split_extend_bed(bed_fn, output_fn, contig_set=None):
    expand_region_size = param.no_of_positions
    output_ctg_dict = defaultdict(list)
    unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (bed_fn)))
    for row_id, row in enumerate(unzip_process.stdout):
        if row[0] == '#':
            continue
        columns = row.strip().split()
        ctg_name = columns[0]
        if contig_set and ctg_name not in contig_set:
            continue

        ctg_start, ctg_end = int(columns[1]), int(columns[2])

        if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
            sys.exit(log_error(
                "[ERROR] Invalid BED input at the {}-th row {} {} {}".format(row_id + 1, ctg_name, ctg_start, ctg_end)))
        expand_ctg_start = max(0, ctg_start - expand_region_size)
        expand_ctg_end = max(0, ctg_end + expand_region_size)
        output_ctg_dict[ctg_name].append(
            ' '.join([ctg_name, str(expand_ctg_start), str(expand_ctg_end)]))

    for key, value in output_ctg_dict.items():
        ctg_output_fn = os.path.join(output_fn, key)
        with open(ctg_output_fn, 'w') as output_file:
            output_file.write('\n'.join(value))

    unzip_process.stdout.close()
    unzip_process.wait()

def write_region_bed(region):

    try:
        ctg_name, start_end = region.split(':')
        ctg_start, ctg_end = int(start_end.split('-')[0]) - 1, int(start_end.split('-')[1]) - 1  # bed format
    except:
        sys.exit("[ERROR] Please use the correct format for --region: ctg_name:start-end, your input is {}".format(
            region))
    if ctg_end < ctg_start or ctg_start < 0 or ctg_end < 0:
        sys.exit("[ERROR] Invalid region input: {}".format(region))

    output_bed_path = os.path.join(args.output_dir, 'tmp', 'region.bed')
    with open(output_bed_path, 'w') as f:
        f.write('\t'.join([ctg_name, str(ctg_start), str(ctg_end)]) + '\n')
    return output_bed_path

def check_contigs_intersection(args, fai_fn):

    MIN_CHUNK_LENGTH = 200000
    MAX_CHUNK_LENGTH = 20000000
    is_include_all_contigs = args.include_all_ctgs
    is_bed_file_provided = args.bed_fn is not None or args.region is not None
    is_known_vcf_file_provided = args.truth_vcf_fn is not None
    is_ctg_name_list_provided = args.ctg_name is not None

    if args.region is not None:
        args.bed_fn = write_region_bed(args.region)

    split_bed_path = os.path.join(args.output_dir, 'tmp', 'split_beds')
    tree = bed_tree_from(bed_file_path=args.bed_fn)
    know_vcf_contig_set = split_extend_vcf(truth_vcf_fn=args.truth_vcf_fn, output_fn=split_bed_path) if is_known_vcf_file_provided else set()
    contig_set = set(args.ctg_name.split(',')) if is_ctg_name_list_provided else set()

    if not args.include_all_ctgs:
        logging("[INFO] --include_all_ctgs not enabled, use chr{1..22,X,Y} and {1..22,X,Y} by default")
    else:
        logging("[INFO] --include_all_ctgs enabled")

    if is_ctg_name_list_provided and is_bed_file_provided:
        logging(log_warning("[WARNING] both --ctg_name and --bed_fn provided, will only proceed with the contigs appeared in both"))

    if is_ctg_name_list_provided and is_known_vcf_file_provided:
        logging(log_warning("[WARNING] both --ctg_name and --truth_vcf_fn provided, will only proceed with the contigs appeared in both"))

    if is_ctg_name_list_provided:
        contig_set = contig_set.intersection(
            set(tree.keys())) if is_bed_file_provided else contig_set
        contig_set = contig_set.intersection(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set
    else:
        contig_set = contig_set.union(
            set(tree.keys())) if is_bed_file_provided else contig_set

        contig_set = contig_set.union(
            know_vcf_contig_set) if is_known_vcf_file_provided else contig_set

    # if each split region is too small(long) for given default chunk num, will increase(decrease) the total chunk num
    default_chunk_num = 0
    DEFAULT_CHUNK_SIZE = args.chunk_size
    contig_length_list = []
    contig_chunk_num = {}

    with open(fai_fn, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name, contig_length = columns[0], int(columns[1])
            if not is_include_all_contigs and (
            not (is_bed_file_provided or is_ctg_name_list_provided or is_known_vcf_file_provided)) and str(
                    contig_name) not in major_contigs:
                continue

            if is_bed_file_provided and contig_name not in tree:
                continue
            if is_ctg_name_list_provided and contig_name not in contig_set:
                continue
            if is_known_vcf_file_provided and contig_name not in contig_set:
                continue

            contig_set.add(contig_name)
            contig_length_list.append(contig_length)
            chunk_num = int(
                contig_length / float(DEFAULT_CHUNK_SIZE)) + 1 if contig_length % DEFAULT_CHUNK_SIZE else int(
                contig_length / float(DEFAULT_CHUNK_SIZE))
            contig_chunk_num[contig_name] = max(chunk_num, 1)

    if default_chunk_num > 0:
        min_chunk_length = min(contig_length_list) / float(default_chunk_num)
        max_chunk_length = max(contig_length_list) / float(default_chunk_num)

    contigs_order = major_contigs_order + list(contig_set)

    sorted_contig_list = sorted(list(contig_set), key=lambda x: contigs_order.index(x))

    found_contig = True
    if not len(contig_set):
        if is_bed_file_provided:
            all_contig_in_bed = ' '.join(list(tree.keys()))
            logging(log_warning("[WARNING] No contig in --bed_fn was found in the reference, contigs in BED {}: {}".format(args.bed_fn, all_contig_in_bed)))
        if is_known_vcf_file_provided:
            all_contig_in_vcf = ' '.join(list(know_vcf_contig_set))
            logging(log_warning("[WARNING] No contig in --truth_vcf_fn was found in the reference, contigs in VCF {}: {}".format(args.truth_vcf_fn, all_contig_in_vcf)))
        if is_ctg_name_list_provided:
            all_contig_in_ctg_name = ' '.join(args.ctg_name.split(','))
            logging(log_warning("[WARNING] No contig in --ctg_name was found in the reference, contigs in contigs list: {}".format(all_contig_in_ctg_name)))
        found_contig = False
    else:
        for c in sorted_contig_list:
            if c not in contig_chunk_num:
                logging(log_warning(("[WARNING] Contig {} given but not found in the reference".format(c))))

        # check contig in bam have support reads
        sorted_contig_list, found_contig = check_contig_in_bam(bam_fn=args.bam_fn, sorted_contig_list=sorted_contig_list,
                                                               samtools=args.samtools)

    if not found_contig:
        log_warning("[WARNING] Exit calling because no contig was found in BAM!")
        sys.exit(0)
    logging('[INFO] Execute RU in contigs: {}'.format(' '.join(sorted_contig_list)))
    logging('[INFO] Number of chunks for each contig: {}'.format(
        ' '.join([str(contig_chunk_num[c]) for c in sorted_contig_list])))

    if default_chunk_num > 0 and max_chunk_length > MAX_CHUNK_LENGTH:
        logging(log_warning(
            '[WARNING] The maximum chunk size set {} is larger than the suggested maximum chunk size {}, consider setting a larger --chunk_num= instead for better parallelism.'.format(
                min_chunk_length, MAX_CHUNK_LENGTH)))

    elif default_chunk_num > 0 and min_chunk_length < MIN_CHUNK_LENGTH:
        logging(log_warning(
            '[WARNING] The minimum chunk size set {} is smaller than the suggested  minimum chunk size {}, consider setting a smaller --chunk_num= instead.'.format(
                min_chunk_length, MIN_CHUNK_LENGTH)))

    if default_chunk_num == 0 and max(contig_length_list) < DEFAULT_CHUNK_SIZE / 5:
        logging(log_warning(
            '[WARNING] The length of the longest contig {} is more than five times smaller smaller than the default chunk size {}, consider setting a smaller --chunk_size= instead for better parallelism.'.format(
                max(contig_length_list), DEFAULT_CHUNK_SIZE)))

    contig_path = os.path.join(args.output_dir, 'tmp', 'CONTIGS')
    with open(contig_path, 'w') as output_file:
        output_file.write('\n'.join(sorted_contig_list))

    chunk_list = []
    chunk_list_path = os.path.join(args.output_dir, 'tmp', 'CHUNK_LIST')
    with open(chunk_list_path, 'w') as output_file:
        for contig_name in sorted_contig_list:
            chunk_num = contig_chunk_num[contig_name] if args.chunk_num is None else args.chunk_num
            for chunk_id in range(1, chunk_num + 1):
                output_file.write(contig_name + ' ' + str(chunk_id) + ' ' + str(chunk_num) + '\n')
                chunk_list.append((contig_name, chunk_id, chunk_num))
    args.chunk_list = chunk_list

    return args


def check_args(args):

    if args.conda_prefix is None:
        if 'CONDA_PREFIX' in os.environ:
            args.conda_prefix = os.environ['CONDA_PREFIX']
        else:
            try:
                python_path = subprocess.run('which python', stdout=subprocess.PIPE, shell=True).stdout.decode().rstrip()
                args.conda_prefix = os.path.dirname(os.path.dirname(python_path))
            except:
                sys.exit(log_error("[ERROR] Conda prefix not found, please activate a correct conda environment."))

    args.bam_fn = file_path_from(file_name=args.bam_fn, exit_on_not_found=True)
    bai_fn = file_path_from(file_name=args.bam_fn, suffix=".bai", exit_on_not_found=False, sep='.')
    crai_fn = file_path_from(file_name=args.bam_fn, suffix=".crai", exit_on_not_found=False, sep='.')
    csi_fn = file_path_from(file_name=args.bam_fn, suffix=".csi", exit_on_not_found=False, sep='.')

    args.ref_fn = file_path_from(file_name=args.ref_fn, exit_on_not_found=True)
    fai_fn = file_path_from(file_name=args.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')

    # args.bed_fn = file_path_from(file_name=args.bed_fn, exit_on_not_found=True)
    args.bed_fn = file_path_from(file_name=args.bed_fn, exit_on_not_found=False)

    args.truth_vcf_fn = file_path_from(file_name=args.truth_vcf_fn, exit_on_not_found=True)

    if bai_fn is None and crai_fn is None and csi_fn is None:
        sys.exit(log_error("[ERROR] BAM index file {} or {} not found. Please run `samtools index $BAM` first.".format(args.bam_fn + '.bai',
                                                                                      args.bam_fn + '.crai')))

    if args.min_af is None:
        args.min_af = param.min_af

    if args.min_coverage is None:
        args.min_coverage = 2

    if args.chunk_size is None:
        args.chunk_size = 5000000

    if args.platform not in {'ont', 'hifi', 'ilmn'}:
        logging(log_error('[ERROR] Invalid platform input, optional: {ont, hifi, ilmn}'))

    if args.skip_steps is not None:
        check_skip_steps_legal(args)

    legal_range_from(param_name="threads", x=args.threads, min_num=1, exit_out_of_range=True)
    legal_range_from(param_name="min_coverage", x=args.min_coverage, min_num=0, exit_out_of_range=True)
    legal_range_from(param_name="min_af", x=args.min_af, min_num=0, max_num=1, exit_out_of_range=True)
    legal_range_from(param_name="chunk_size", x=args.chunk_size, min_num=0, exit_out_of_range=True)

    args.output_path = create_output_folder(args)
    check_tools_version(args=args)
    args = check_contigs_intersection(args=args, fai_fn=fai_fn)

    return args


def print_args(args):
    version = '0.0.1'
    logging("")
    logging("[INFO] RU VERSION: {}".format(version))
    logging("[INFO] BAM FILE PATH: {}".format(args.bam_fn))
    logging("[INFO] REFERENCE FILE PATH: {}".format(args.ref_fn))
    logging("[INFO] TRUTH VCF FILE PATH: {}".format(args.truth_vcf_fn))
    logging("[INFO] PLATFORM: {}".format(args.platform))
    logging("[INFO] THREADS: {}".format(args.threads))
    logging("[INFO] OUTPUT FOLDER: {}".format(args.output_dir))
    logging("[INFO] OUTPUT VCF PATH: {}".format(os.path.join(args.output_dir, 'unified.vcf.gz')))
    logging("[INFO] BED FILE PATH: {}".format(args.bed_fn))
    logging("[INFO] REGION FOR EXECUTING: {}".format(args.region))
    logging("[INFO] CONTIGS FOR EXECUTING: {}".format(args.ctg_name))
    logging("[INFO] CONDA BINARY PREFIX: {}".format(args.conda_prefix))
    logging("[INFO] SAMTOOLS BINARY PATH: {}".format(args.samtools))
    logging("[INFO] WHATSHAP BINARY PATH: {}".format(args.whatshap))
    logging("[INFO] PYTHON BINARY PATH: {}".format(args.python))
    logging("[INFO] PYPY BINARY PATH: {}".format(args.pypy))
    logging("[INFO] PARALLEL BINARY PATH: {}".format(args.parallel))
    logging("[INFO] CHUNK SIZE: {}".format(args.chunk_size))
    logging("[INFO] MINIMUM AF: {}".format(args.min_af))
    logging("[INFO] DISABLE PHASING: {}".format(args.disable_phasing))
    logging("[INFO] ENABLE DRY RUN: {}".format(args.dry_run))
    logging("[INFO] ENABLE INCLUDING ALL CTGS FOR EXECUTING: {}".format(args.include_all_ctgs))
    logging("[INFO] ENABLE REMOVING INTERMEDIATE FILES: {}".format(args.remove_intermediate_dir))
    logging("")

    if args.cmdline is not None and args.cmdline != "":
        with open(args.output_dir + '/tmp/CMD', 'w') as f:
            f.write(args.cmdline + '\n')

    return args


def print_command_line(args):

    try:

        cmdline = os.path.realpath(__file__)
        cmdline += ' --bam_fn {} '.format(args.bam_fn)
        cmdline += '--ref_fn {} '.format(args.ref_fn)
        cmdline += '--truth_vcf_fn {} '.format(args.truth_vcf_fn)
        cmdline += '--threads {} '.format(args.threads)
        cmdline += '--platform {} '.format(args.platform)
        cmdline += '--output_dir {} '.format(args.output_dir)
        cmdline += '--bed_fn {} '.format(args.bed_fn) if args.bed_fn is not None else ""
        cmdline += '--ctg_name {} '.format(args.ctg_name) if args.ctg_name is not None else ""
        cmdline += '--region {} '.format(args.region) if args.region is not None else ""
        cmdline += '--min_af {} '.format(args.snv_min_af) if args.min_af is not None else ""
        cmdline += '--min_coverage {} '.format(args.min_coverage) if args.min_coverage is not None else ""
        cmdline += '--chunk_size {} '.format(args.chunk_size) if args.chunk_size is not None else ""
        cmdline += '--sample_name {} '.format(args.sample_name) if args.sample_name != "SAMPLE" else ""
        cmdline += '--output_prefix {} '.format(args.output_prefix) if args.output_prefix != "output" else ""
        cmdline += '--remove_intermediate_dir ' if args.remove_intermediate_dir else ""
        cmdline += '--include_all_ctgs ' if args.include_all_ctgs else ""
        cmdline += '--python {} '.format(args.python) if args.python != "python3" else ""
        cmdline += '--pypy {} '.format(args.pypy) if args.pypy != "pypy3" else ""
        cmdline += '--samtools {} '.format(args.samtools) if args.samtools != "samtools" else ""
        cmdline += '--whatshap {} '.format(args.whatshap) if args.whatshap != "whatshap" else ""
        cmdline += '--parallel {} '.format(args.parallel) if args.parallel != "parallel" else ""
        cmdline += '--disable_phasing ' if args.disable_phasing else ""
        cmdline += '--skip_steps {} '.format(args.skip_steps) if args.skip_steps is not None else ""
        cmdline += '--conda_prefix {} '.format(args.conda_prefix) if args.conda_prefix is not None else ""
        args.cmdline = cmdline
    except:
        return args
    logging("[COMMAND] " + cmdline + '\n')
    return args


def ru(args):

    step = 1
    echo_list = []
    commands_list = []

    try:
        rc = subprocess.check_call('time', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time = 'time '
    except subprocess.CalledProcessError as e:
        time = ''

    # RU
    # STEP 1: Whatshap UnPhasing
    echo_list.append("[INFO] STEP {}: Whatshap UnPhasing".format(step))
    step += 1
    wu_command = args.whatshap + ' unphase {} > {}'.format(args.truth_vcf_fn, args.output_dir + '/tmp/INPUT.vcf.gz')
    commands_list += [wu_command]

    # STEP 2: Whatshap Phasing
    echo_list.append("[INFO] STEP {}: Whatshap Phasing INPUT VCF file".format(step))
    step += 1
    wp_command = '( ' + time + args.parallel
    wp_command += ' --joblog ' + args.output_dir + '/logs/parallel_1_phase.log'
    wp_command += ' -j ' + str(args.threads)
    wp_command += ' ' + args.whatshap + ' phase'
    wp_command += ' --output ' + args.output_dir + '/tmp/phased_vcf/phased_{1}.vcf.gz'
    wp_command += ' --reference ' + args.ref_fn
    wp_command += ' --chromosome {1}'
    wp_command += ' --ignore-read-groups'
    wp_command += ' --distrust-genotypes'
    wp_command += ' ' + args.output_dir + '/tmp/INPUT.vcf.gz'
    wp_command += ' ' + args.bam_fn
    wp_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    wp_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/1_WP.log'
    commands_list += [wp_command]

    # STEP 3: Index VCF
    echo_list.append("[INFO] STEP {}: Index VCF file".format(step))
    step += 1
    tabix_input_command = args.parallel + ' -j ' + str(args.threads)
    tabix_input_command += ' tabix' + ' -f -p vcf'
    tabix_input_command += ' ' + args.output_dir + '/tmp/phased_vcf/phased_{1}.vcf.gz'
    tabix_input_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    commands_list.append(tabix_input_command)

    # STEP 4: Whatshap Haplotagging
    echo_list.append("[INFO] STEP {}: Whatshap Haplotagging BAM file".format(step))
    step += 1
    wh_command = '( ' + time + args.parallel
    wh_command += ' --joblog ' + args.output_dir + '/logs/parallel_2_haplotag.log'
    wh_command += ' -j ' + str(args.threads)
    wh_command += ' ' + args.whatshap + ' haplotag'
    wh_command += ' --output ' + args.output_dir + '/tmp/phased_bam/{2}_{1}.bam'
    wh_command += ' --reference ' + args.ref_fn
    wh_command += ' --regions {1}'
    wh_command += ' --ignore-read-groups'
    wh_command += ' ' + args.output_dir + '/tmp/phased_vcf/phased_{1}.vcf.gz'
    wh_command += ' ' + args.bam_fn
    wh_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    wh_command += ' ::: ' + args.sample_name
    wh_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/2_WH.log'

    index_command = args.parallel + ' -j ' + str(args.threads)
    index_command += ' ' + args.samtools + ' index '
    index_command += ' -@' + str(args.threads)
    index_command += ' ' + args.output_dir + '/tmp/phased_bam/{2}_{1}.bam'
    index_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    index_command += ' ::: ' + args.sample_name

    commands_list.append(wh_command + ' && ' + index_command)

    # STEP 5: Split bed file regions according to the contig name and extend bed region
    echo_list.append("[INFO] STEP {}: Split bed file regions according to the contig name and extend bed region".format(step))
    step += 1
    seb_command = '( ' + time + args.parallel
    seb_command += ' --joblog ' + args.output_dir + '/logs/parallel_3_split_extend_bed.log'
    seb_command += ' -j ' + str(args.threads)
    seb_command += ' ' + args.pypy + ' ' + main_entry + ' SplitExtendBed'
    seb_command += ' --bed_fn {}'.format(args.bed_fn) if args.bed_fn is not None else ''
    seb_command += ' --ctgName {1}'
    seb_command += ' --output_fn ' + args.output_dir + '/tmp/split_beds/{1}'
    seb_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    seb_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/3_SEB.log'
    commands_list.append(seb_command)

    # STEP 6: Get true variant label information from VCF file
    echo_list.append("[INFO] STEP {}: Get true variant label information from VCF file".format(step))
    step += 1
    gtv_command = '( ' + time + args.parallel
    gtv_command += ' --joblog ' + args.output_dir + '/logs/parallel_4_get_true_variant.log'
    gtv_command += ' -j ' + str(args.threads)
    gtv_command += ' ' + args.pypy + ' ' + main_entry + ' GetTruth'
    gtv_command += ' --vcf_fn ' + args.truth_vcf_fn
    gtv_command += ' --ctgName {1}'
    gtv_command += ' --var_fn ' + args.output_dir + '/tmp/var/var_{1}'
    gtv_command += ' :::: ' + args.output_dir + '/tmp/CONTIGS'
    gtv_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/4_GTV.log'
    commands_list.append(gtv_command)

    # STEP 7: Extract Variant Candidates from BAM file
    echo_list.append("[INFO] STEP {}: Extract variant candidates from BAM file".format(step))
    step += 1
    ec_command = '( ' + time + args.parallel
    ec_command += ' --joblog ' + args.output_dir + '/logs/parallel_5_extract_candidates.log'
    ec_command += ' -C " " -j ' + str(args.threads)
    ec_command += ' ' + args.pypy + ' ' + main_entry + ' ExtractCandidates'
    ec_command += ' --bam_fn ' + args.output_dir + '/tmp/phased_bam/{4}_{1}.bam'
    ec_command += ' --ref_fn ' + args.ref_fn
    ec_command += ' --ctgName {1}'
    ec_command += ' --samtools ' + args.samtools
    ec_command += ' --min_af ' + str(args.min_af)
    ec_command += ' --extend_bed ' + args.output_dir + '/tmp/split_beds/{1}'
    ec_command += ' --chunk_id {2}'
    ec_command += ' --chunk_num {3}'
    ec_command += ' --phasing_info_in_bam'
    ec_command += ' --test_pos False'
    ec_command += ' --unify_repre_fn ' + args.output_dir + '/tmp/candidates/{4}_{2}_{1}'
    ec_command += ' :::: ' + os.path.join(args.output_dir, 'tmp', 'CHUNK_LIST')
    ec_command += ' ::: ' + args.sample_name
    ec_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/5_EC.log'
    commands_list.append(ec_command)

    # STEP 8: Label Variant
    echo_list.append("[INFO] STEP {}: Label variant".format(step))
    step += 1
    lv_command = '( ' + time + args.parallel
    lv_command += ' --joblog ' + args.output_dir + '/logs/parallel_6_unify_representation.log'
    lv_command += ' -C " " -j ' + str(args.threads)
    lv_command += ' ' + args.pypy + ' ' + main_entry + ' UnifyRepresentation'
    lv_command += ' --vcf_fn ' + args.output_dir + '/tmp/var/var_{1}'
    lv_command += ' --candidates_fn ' + args.output_dir + '/tmp/candidates/{4}_{2}_{1}'
    lv_command += ' --ref_fn ' + args.ref_fn
    lv_command += ' --output_vcf ' + args.output_dir + '/tmp/vcf_output/vcf_{1}_{2}'
    lv_command += ' --platform ' + args.platform
    lv_command += ' --ctgName {1}'
    lv_command += ' --bed_fn ' + args.output_dir + '/tmp/split_beds/{1}'
    lv_command += ' --chunk_id {2}'
    lv_command += ' --chunk_num {3}'
    lv_command += ' --test_pos False'
    lv_command += ' --sampleName {4}'
    lv_command += ' --subsample_ratio 1000'
    lv_command += ' :::: ' + os.path.join(args.output_dir, 'tmp', 'CHUNK_LIST')
    lv_command += ' ::: ' + args.sample_name
    lv_command += ' ) 2>&1 | tee ' + args.output_dir + '/logs/6_RU.log'
    commands_list.append(lv_command)

    # STEP 9: Sort VCF
    echo_list.append("[INFO] STEP {}: Sort VCF file".format(step))
    step += 1
    sv_command = 'cat ' + args.output_dir + '/tmp/vcf_output/vcf_*'
    sv_command += ' | '
    sv_command += args.pypy + ' ' + main_entry + ' SortVcf'
    sv_command += ' --output_fn ' + args.output_dir + '/unified.vcf'
    sv_command += ' && '
    sv_command += 'bgzip' + ' -f ' + args.output_dir + '/unified.vcf'
    sv_command += ' && '
    sv_command += 'tabix' + ' -f -p vcf ' + args.output_dir + '/unified.vcf.gz'
    commands_list.append(sv_command)

    # Execute commands step by step
    skip_steps = args.skip_steps.rstrip().split(',') if args.skip_steps else None
    stdout = sys.stdout if args.tee is None else args.tee.stdin
    for i, (command, echo) in enumerate(zip(commands_list, echo_list)):

        logging(echo)
        logging("[INFO] RUN THE FOLLOWING COMMAND:")
        logging(command)
        logging("")
        if not args.dry_run:
            if skip_steps is not None and str(i+1) in skip_steps:
                logging("[INFO] --skip_steps is enabled, skip running step {}.".format(i+1))
                logging("")
                continue
            try:
                return_code = subprocess.check_call(command, shell=True, stdout=stdout)
            except subprocess.CalledProcessError as e:
                sys.stderr.write("ERROR in STEP {}, THE FOLLOWING COMMAND FAILED: {}\n".format(i+1, command))
                exit(1)
        logging("")

    if args.remove_intermediate_dir:
        logging("[INFO] Removing intermediate files in {}/tmp ...".format(args.output_dir))
        subprocess.run('rm -rf {}/tmp'.format(args.output_dir), shell=True)


def ru_parser():

    parser = argparse.ArgumentParser(
        description="Run Clair3 Representation Unification. Example run: ru -b BAM_FN -r REF_FN --truth_vcf_fn TRUTH_VCF_FN -o OUTPUT_DIR -t THREADS -p PLATFORM")

    # print version
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format('0.0.1'))

    required_params = parser.add_argument_group('Required parameters')
    required_params.add_argument(
        '-b',
        "--bam_fn",
        type=str,
        required=True,
        default=None,
        help="BAM file input. The input file must be samtools indexed."
    )

    required_params.add_argument(
        "-r",
        "--ref_fn",
        type=str,
        required=True,
        default=None,
        help="FASTA reference file input. The input file must be samtools indexed."
    )

    required_params.add_argument(
        "--truth_vcf_fn",
        type=str,
        required=True,
        default=None,
        help="Truth VCF file input."
    )

    required_params.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        default=None,
        help="Output directory."
    )

    required_params.add_argument(
        "-t",
        "--threads",
        required=True,
        type=int,
        default=None,
        help="Max #threads to be used."
    )

    required_params.add_argument(
        "-p",
        "--platform",
        required=True,
        type=str,
        default=None,
        help="Select the sequencing platform of the input. Possible options: {ont, hifi, ilmn}."
    )

    optional_params = parser.add_argument_group('Optional parameters')
    optional_params.add_argument(
        "-c",
        "--ctg_name",
        type=str,
        default=None,
        help="The name of the contigs to be processed. Split by ',' for multiple contigs. Default: all contigs will be processed."
    )

    optional_params.add_argument(
        "--bed_fn",
        type=str,
        default=None,
        help="Path to a BED file. Execute RU only in the provided BED regions."
    )

    optional_params.add_argument(
        "--region",
        type=str,
        default=None,
        help="A region to be processed. Format: `ctg_name:start-end` (start is 1-based)."
    )

    optional_params.add_argument(
        "--min_af",
        type=float,
        default=None,
        help="Minimal AF required for a variant to be called. Default: 0.08."
    )

    optional_params.add_argument(
        "--min_coverage",
        type=int,
        default=None,
        help="Minimal coverage required for a variant to be called. Default: 4."
    )

    optional_params.add_argument(
        "-s",
        "--sample_name",
        type=str,
        default="SAMPLE",
        help="Define the sample name to be shown in the VCF file. Default: SAMPLE."
    )

    optional_params.add_argument(
        "--output_prefix",
        type=str,
        default="output",
        help="Prefix for output VCF filename. Default: output."
    )

    optional_params.add_argument(
        "--remove_intermediate_dir",
        action='store_true',
        help="Remove intermediate directory before finishing to save disk space."
    )

    optional_params.add_argument(
        "--include_all_ctgs",
        action='store_true',
        help="Execute RU on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}."
    )

    optional_params.add_argument(
        '-d',
        "--dry_run",
        action='store_true',
        help="Print the commands that will be ran."
    )

    optional_params.add_argument(
        "--python",
        type=str,
        default="python3",
        help="Absolute path of python, python3 >= 3.9 is required."
    )

    optional_params.add_argument(
        "--pypy",
        type=str,
        default="pypy3",
        help="Absolute path of pypy3, pypy3 >= 3.6 is required."
    )

    optional_params.add_argument(
        "--samtools",
        type=str,
        default="samtools",
        help="Absolute path of samtools, samtools version >= 1.10 is required."
    )

    optional_params.add_argument(
        "--whatshap",
        type=str,
        default="whatshap",
        help="Absolute path of whatshap, whatshap >= 1.0 is required."
    )

    optional_params.add_argument(
        "--parallel",
        type=str,
        default="parallel",
        help="Absolute path of parallel, parallel >= 20191122 is required."
    )

    optional_params.add_argument(
        "--disable_phasing",
        action='store_true',
        help="Disable phasing with whatshap."
    )

    optional_params.add_argument(
        "--chunk_size",
        type=int,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--chunk_num",
        type=int,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--output_path",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--skip_steps",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--tee",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--conda_prefix",
        type=str,
        default=None,
        help=SUPPRESS
    )

    optional_params.add_argument(
        "--cmdline",
        type=str,
        default=None,
        help=SUPPRESS
    )

    return parser

def main():
    """
    Main interface for Clair3 RU.
    """

    global args

    call_start_time = time()

    parser = ru_parser()
    args = parser.parse_args()

    args.output_dir = folder_path_from(args.output_dir, create_not_found=True)
    tee_logger = os.path.join(args.output_dir, 'ru.log' if not args.dry_run else "ru_dry_run.log")
    if os.path.exists(tee_logger):
        subprocess.run("mv {} {}".format(tee_logger, tee_logger + '.bak'), shell=True)
    try:
        args.tee = subprocess.Popen(['tee', tee_logger], stdin=subprocess.PIPE, bufsize=0)
    except:
        logging(log_warning("[WARNING] `tee` not found, disable `tee` logging!"))
        args.tee = None

    logging("")

    args = print_command_line(args)
    args = check_args(args)
    args = print_args(args)
    ru(args)

    runtime = time() - call_start_time
    logging("[INFO] Total time elapsed: %im%.2fs\n" % (int(runtime/60), int(runtime % 60)))
    logging("[INFO] Finish Representation Unification, output vcf file path: {}/unified.vcf.gz\n".format(args.output_dir))

    if args.tee is not None:
        args.tee.stdin.close()

if __name__ == '__main__':
    main()
