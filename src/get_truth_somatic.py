import sys
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from
import os
import subprocess
class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()

def OutputVariant(args):
    var_fn = args.var_fn
    vcf_fn = args.vcf_fn
    truth_vcf_fn = args.truth_vcf_fn
    ctg_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd

    truth_vcf_set = set()
    variant_set = set()
    compress = False
    if args.truth_vcf_fn is not None:
        truth_vcf_set = set(vcf_candidates_from(vcf_fn=truth_vcf_fn, contig_name=ctg_name))
    if args.var_fn != "PIPE":
        if compress:
            var_fpo = open(var_fn, "wb")
            var_fp = subprocess_popen(shlex.split("gzip -c"), stdin=PIPE, stdout=var_fpo)
        else:
            var_fpo = open(var_fn, "w")
            var_fp = var_fpo
    else:
        var_fp = TruthStdout(sys.stdout)

    is_ctg_region_provided = ctg_start is not None and ctg_end is not None

    vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))

    if args.output_bed_fn is not None:
        output_dir = os.path.dirname(args.output_bed_fn)
        if not os.path.exists(output_dir):
            subprocess.run("mkdir -p {}".format(output_dir), shell=True)
        output_bed_file = open(args.output_bed_fn, 'w')

    for row in vcf_fp.stdout:
        columns = row.strip().split()
        if columns[0][0] == "#":
            continue

        # position in vcf is 1-based
        chromosome, position = columns[0], columns[1]
        if chromosome != ctg_name:
            continue
        if is_ctg_region_provided and not (ctg_start <= int(position) <= ctg_end):
            continue
        reference, alternate, last_column = columns[3], columns[4], columns[-1]


        af_index = columns[8].split(':').index('AF')
        af = columns[9].split(':')[af_index]
        if float(af) > args.max_af:
            continue
        genotype_1, genotype_2 = '0', '1'

        variant_set.add(int(position))
        if not compress:
            var_fp.write(" ".join((chromosome, position, reference, alternate, genotype_1, genotype_2)) + "\n")
        else:
            var_fp.stdin.write(" ".join((chromosome, position, reference, alternate, genotype_1, genotype_2)) + "\n")

        if args.output_bed_fn is not None:
            output_bed_file.write("{}\t{}\t{}\n".format(chromosome, int(position)-50, int(position)+50))

    for position in truth_vcf_set:
        if position not in variant_set:
            # miss variant set used in Tensor2Bin
            if not compress:
                var_fp.write(" ".join((chromosome, str(position), "None", "None", "-1", "-1")) + "\n")
            else:
                var_fp.stdin.write(" ".join((chromosome, str(position), "None", "None", "-1", "-1")) + "\n")


    vcf_fp.stdout.close()
    vcf_fp.wait()

    if args.var_fn != "PIPE":
        if compress:
            var_fp.stdin.close()
            var_fp.wait()
            var_fpo.close()
        else:
            var_fp.close()

def main():
    parser = ArgumentParser(description="Extract variant type and allele from a truth dataset")

    parser.add_argument('--vcf_fn', type=str, default="input.vcf", required=True,
                        help="Truth VCF file input, required")

    parser.add_argument('--var_fn', type=str, default="PIPE",
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--output_bed_fn', type=str, default=None,
                        help="Truth variants output, use PIPE for standard output, default: %(default)s")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")

    parser.add_argument('--max_af', type=float, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed")


    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Truth VCF file input, only used when vcf_fn is unified vcf. Marked truth variants not in unified as missing")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    OutputVariant(args)


if __name__ == "__main__":
    main()
