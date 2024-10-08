import sys
import shlex
import os

from textwrap import dedent
from subprocess import PIPE, run
from argparse import ArgumentParser
from collections import defaultdict
from shared.utils import subprocess_popen, Position as Position, file_path_from


class TruthStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()

class VcfWriter(object):
    def __init__(self, vcf_fn, ctg_name=None, ref_fn=None, sample_name="SAMPLE", write_header=True, show_ref_calls=False):
        self.vcf_fn = vcf_fn
        self.show_ref_calls = show_ref_calls
        #make directory if not exist
        vcf_folder = os.path.dirname(self.vcf_fn)
        if not os.path.exists(vcf_folder):
            print("[INFO] Output VCF folder {} not found, create it".format(vcf_folder))
            return_code = run("mkdir -p {}".format(vcf_folder), shell=True)

        self.vcf_writer = open(self.vcf_fn, 'w')
        self.ref_fn = ref_fn
        self.ctg_name = ctg_name
        self.sample_name = sample_name
        if write_header:
            self.write_header(ref_fn=ref_fn)

    def close(self):
        try:
            self.vcf_writer.close()
        except:
            pass

    def write_header(self, ctg_name=None, ref_fn=None):
        header = dedent("""\
                    ##fileformat=VCFv4.2
                    ##FILTER=<ID=PASS,Description="All filters passed">
                    ##FILTER=<ID=LowQual,Description="Low quality variant">
                    ##FILTER=<ID=RefCall,Description="Reference call">
                    ##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
                    ##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
                    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Overall Read Depth(Normal+Tumor)">
                    ##FORMAT=<ID=NDP,Number=1,Type=Integer,Description="Normal Read Depth">
                    ##FORMAT=<ID=TDP,Number=1,Type=Integer,Description="Tumor Read Depth">
                    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
                    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
                    ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range of [0,1]">
                    ##FORMAT=<ID=AU,Number=1,Type=Integer,Description="Number of 'A' alleles in Tumor BAM">
                    ##FORMAT=<ID=CU,Number=1,Type=Integer,Description="Number of 'C' alleles in Tumor BAM">
                    ##FORMAT=<ID=GU,Number=1,Type=Integer,Description="Number of 'G' alleles in Tumor BAM">
                    ##FORMAT=<ID=TU,Number=1,Type=Integer,Description="Number of 'T' alleles in Tumor BAM">
                    """
               )

        if self.ref_fn is not None:
            reference_index_file_path = file_path_from(self.ref_fn, suffix=".fai", exit_on_not_found=True, sep='.')
            with open(reference_index_file_path, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    if self.ctg_name is not None and contig_name != self.ctg_name:
                        continue
                    header += "##contig=<ID=%s,length=%s>\n" % (contig_name, contig_size)

        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (self.sample_name)

        self.vcf_writer.write(header)

    def write_row(self, POS=None, REF=None, ALT=None, QUAL=0, GT='0/0', DP=0, AF=0, CHROM=None, GQ=None, ID='.', FILTER=".", INFO='.', NAF=None, TAF=None, VT=None,
                  NDP=None, TDP=None, AU=None, CU=None, GU=None, TU=None, row_str=None):
        if row_str is not None:
            self.vcf_writer.write(row_str)
            return
        GQ = GQ if GQ else QUAL
        CHROM = CHROM if CHROM else self.ctg_name
        if not self.show_ref_calls and (GT == "0/0" or GT == "./."):
            return
        FORMAT = "GT:GQ:DP:AF"
        FORMAT_V = "%s:%d:%d:%.4f" % (GT, GQ, DP, AF)
        basic_vcf_format = "%s\t%d\t%s\t%s\t%s\t%.4f\t%s\t%s" % (
            CHROM,
            int(POS),
            ID,
            REF,
            ALT,
            QUAL,
            FILTER,
            INFO
            )
        if NAF is not None:
            FORMAT += ":NAF"
            FORMAT_V += ":%.4f" % (NAF)
        if TAF is not None:
            FORMAT += ":TAF"
            FORMAT_V += ":%.4f" % (TAF)
        if NDP is not None:
            FORMAT += ":NDP"
            FORMAT_V += ":%d" % (NDP)
        if TDP is not None:
            FORMAT += ":TDP"
            FORMAT_V += ":%d" % (TDP)
        if AU is not None and CU is not None and GU is not None and TU is not None:
            FORMAT += ":AU:CU:GU:TU"
            FORMAT_V += ":%d:%d:%d:%d" % (AU, CU, GU, TU)
        if VT is not None:
            FORMAT += ":VT"
            FORMAT_V += ":%s" % (VT)
        vcf_format = '\t'.join([basic_vcf_format, FORMAT, FORMAT_V]) + "\n"

        self.vcf_writer.write(vcf_format)

class VcfReader(object):
    def __init__(self, vcf_fn,
                 ctg_name=None,
                 ctg_start=None,
                 ctg_end=None,
                 is_var_format=False,
                 is_happy_format=False,
                 is_fp=None,
                 show_ref=True,
                 direct_open=False,
                 keep_row_str=False,
                 skip_genotype=False,
                 filter_tag=None,
                 naf_filter=None,
                 taf_filter=None,
                 save_header=False,
                 min_qual=None,
                 max_qual=None,
                 discard_indel=False,
                 keep_af=False):
        self.vcf_fn = vcf_fn
        self.ctg_name = ctg_name
        self.ctg_start = ctg_start
        self.ctg_end = ctg_end
        self.variant_dict = defaultdict(Position)
        self.is_var_format = is_var_format
        self.is_happy_format = is_happy_format
        self.is_fp = is_fp
        self.show_ref = show_ref
        self.direct_open = direct_open
        self.keep_row_str = keep_row_str
        self.skip_genotype = skip_genotype
        self.filter_tag = filter_tag #PASS;HighConf PASS;MedConf in hcc1395
        self.naf_filter = naf_filter
        self.taf_filter = taf_filter
        self.header = ""
        self.save_header = save_header
        self.discard_indel = discard_indel
        self.min_qual = min_qual
        self.max_qual = max_qual
        self.keep_af = keep_af
    def read_vcf(self):
        is_ctg_region_provided = self.ctg_start is not None and self.ctg_end is not None

        if self.vcf_fn is None or not os.path.exists(self.vcf_fn):
            return

        vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (self.vcf_fn)))
        for row in vcf_fp.stdout:
            columns = row.strip().split()
            if columns[0][0] == "#":
                continue

            # position in vcf is 1-based
            chromosome, position = columns[0], columns[1]
            if chromosome != self.ctg_name:
                continue
            if is_ctg_region_provided and not (self.ctg_start <= int(position) <= self.ctg_end):
                continue
            if self.is_var_format:
                reference, alternate = columns[2], columns[3]
                genotype_1 = int(columns[4])
                genotype_2 = int(columns[5])
            else:
                reference, alternate, last_column = columns[3], columns[4], columns[-1]
            # normal GetTruth
                if self.is_happy_format and self.is_fp:
                    last_column = columns[10]
                if self.is_happy_format and not self.is_fp:
                    last_column = columns[9]
                try:
                    genotype = last_column.split(":")[0].replace("/", "|").replace(".", "0").split("|")
                    if len(genotype) == 1:
                        genotype_1 = genotype[0]
                        genotype_2 = 3
                    else:
                        genotype_1, genotype_2 = genotype

                    # 1000 Genome GetTruth (format problem) (no genotype is given)
                    if int(genotype_1) > int(genotype_2):
                        genotype_1, genotype_2 = genotype_2, genotype_1

                    #remove * to guarentee vcf match
                    if '*' in alternate:
                        alternate = alternate.split(',')
                        if int(genotype_1) + int(genotype_2) != 3 or len(alternate) != 2:
                            print ('error with variant representation')
                            continue
                        alternate = ''.join([alt_base for alt_base in alternate if alt_base != '*'])
                        # * always have a genotype 1/2

                        genotype_1, genotype_2 = '0', '1'
                except:
                    genotype_1, genotype_2 = '0', '0'
            position = int(position)
            self.variant_dict[position] = Position(pos=position,
                                                    ref_base=reference,
                                                   alt_base=alternate,
                                                   genotype1=int(genotype_1),
                                                   genotype2=int(genotype_2))
    def get_alt_info(self, pos, extra_info=""):
        pos = int(pos)
        if pos not in self.variant_dict:
            return ""
        ref_base = self.variant_dict[pos].reference_bases
        alt_base = ','.join(self.variant_dict[pos].alternate_bases)
        gentoype_str = '/'.join([str(g) for g in self.variant_dict[pos].genotype])
        return extra_info + '_' + ref_base + '_' + alt_base + '_' + gentoype_str