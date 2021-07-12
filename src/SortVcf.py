import os
import subprocess

from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict

from shared.utils import log_error, log_warning, file_path_from, output_header
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]


def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



def print_calling_step(output_fn=""):

    merge_output = os.path.join(os.path.dirname(output_fn), 'merge_output.vcf.gz')
    pileup_output = os.path.join(os.path.dirname(output_fn), 'pileup.vcf.gz')

    print (log_warning("[WARNING] Copying pileup.vcf.gz to {}".format(merge_output)))
    subprocess.run('cp {} {}'.format(pileup_output, merge_output), shell=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE)

def sort_vcf_from_stdin(args):
    """
    Sort vcf file according to variants start position and contig name.
    """
    
    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    for row in stdin:
        row_count += 1
        if row[0] == '#':
            if row not in header:
                header.append(row)
            continue
        # use the first vcf header
        columns = row.strip().split(maxsplit=3)
        ctg_name, pos = columns[0], columns[1]
        contig_dict[ctg_name][int(pos)] = row
        no_vcf_output = False
    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, please check the setting"))
    if no_vcf_output:
        print(log_warning("[WARNING] No variant found, please check the setting"))

    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(args.output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])


def sort_vcf_from(args):
    """
    Sort vcf file from providing vcf filename prefix.
    """
    output_fn = args.output_fn
    input_dir = args.input_dir
    vcf_fn_prefix = args.vcf_fn_prefix
    vcf_fn_suffix = args.vcf_fn_suffix
    sample_name= args.sampleName
    ref_fn= args.ref_fn

    if not os.path.exists(input_dir):
        exit(log_error("[ERROR] Input directory: {} not exists!").format(input_dir))
    all_files = os.listdir(input_dir)

    if vcf_fn_prefix is not None:
        all_files = [item for item in all_files if item.startswith(vcf_fn_prefix)]
        if len(all_files) == 0:
            output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
            print (log_warning(
                "[WARNING] No vcf file found with prefix:{}/{}, output empty vcf file".format(input_dir,vcf_fn_prefix)))
            compress_index_vcf(output_fn)
            print_calling_step(output_fn=output_fn)
            return

    if vcf_fn_suffix is not None:
        all_files = [item for item in all_files if item.endswith(vcf_fn_suffix)]
        if len(all_files) == 0:
            output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
            print (log_warning(
                "[WARNING] No vcf file found with suffix:{}/{}, output empty vcf file".format(input_dir,vcf_fn_prefix)))
            compress_index_vcf(output_fn)
            print_calling_step(output_fn=output_fn)
            return

    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    for vcf_fn in all_files:
        for row in open(os.path.join(input_dir, vcf_fn), 'r'):
            row_count += 1
            if row[0] == '#':
                if row not in header:
                    header.append(row)
                continue
            # use the first vcf header
            columns = row.strip().split(maxsplit=3)
            ctg_name, pos = columns[0], columns[1]
            contig_dict[ctg_name][int(pos)] = row
            no_vcf_output = False
    if row_count == 0:
        print (log_warning("[WARNING] No vcf file found, output empty vcf file"))
        output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
        compress_index_vcf(output_fn)
        print_calling_step(output_fn=output_fn)
        return
    if no_vcf_output:
        output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
        print (log_warning("[WARNING] No variant found, output empty vcf file"))
        compress_index_vcf(output_fn)
        print_calling_step(output_fn=output_fn)
        return

    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])

    compress_index_vcf(output_fn)


def main():
    parser = ArgumentParser(description="Sort a VCF file according to contig name and starting position")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Output VCF filename, required")

    parser.add_argument('--input_dir', type=str, default=None,
                        help="Input directory")

    parser.add_argument('--vcf_fn_prefix', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--vcf_fn_suffix', type=str, default='.vcf',
                        help="Input vcf filename suffix")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    args = parser.parse_args()
    if args.input_dir is None:
        sort_vcf_from_stdin(args)
    else:
        sort_vcf_from(args)
if __name__ == "__main__":
    main()
