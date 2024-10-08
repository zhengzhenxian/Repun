import collections
import heapq
import itertools
import json
import shlex
import os

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from sys import stderr
from subprocess import PIPE, Popen
from shared.vcf import VcfReader

import shared.param as param

from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, region_from, reference_sequence_from, Position as Position, output_header, str2bool
from src.utils import file_path_from
extended_window_size = 200
region_size = 50000
reference_region_size = region_size * 2
extend_bp = 100
reference_allele_gap = 0.9

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

class Reference(object):
    """
    Reference region query with given reference start and reference end, we cocnat the reference base with altertive base
    with reference base to generate read-level haplotype.
    """
    def __init__(self, seq, start, reference_sequence, reference_start):
        self.start = start
        self.end = start + len(seq)
        self.seq = seq
        self.reference_sequence = reference_sequence
        self.reference_start = reference_start

    def query(self, start, end):
        return self.seq[start - self.start:end - self.start]

def unique_genotypes_selection(genotype_options):
    """
    Extend true two haplotypes according to the chosen genotype and only save the haplotype set with distinct haplotype,
    if two haplotypes set match in either hap1(1) == hap2(1) or hap1(1) = hap2(1), remove the duplication to reduce
    massive comparison.
    """

    genotypes_list = []
    existed_genotypes_path = set()
    for genotype_pair in itertools.product(*genotype_options):
        g1, g2 = "", ""
        for h1, h2 in genotype_pair:
            g1 += str(h1)
            g2 += str(h2)
        genotype_tupe = (g1, g2)
        genotype_tupe_reverse = (g2, g1)
        if genotype_tupe in existed_genotypes_path or genotype_tupe_reverse in existed_genotypes_path:
            continue
        existed_genotypes_path.add(genotype_tupe)
        existed_genotypes_path.add(genotype_tupe_reverse)
        genotypes_list.append(genotype_pair)
    return genotypes_list


def find_read_support(variants, ref, variant_type, max_calculate_count, read_seqs_counter=None, candidate_genotypes=None, variant_dict=None, read_name_info_dict=None, truths=None, alt_dict=None, no_match_found=False, output_seq=False):
    """
    Find read-level support for each matched haplotype, we only extended the reference sequence with the alternative base,
    and discard low allele frequency systematic error.
    """
    if variant_type == 'candidate':
        all_read_set = set()
        all_read_seq_dict = defaultdict(str)
        for v in variants:
            pos = v.start
            read_name_set = alt_dict[pos].read_name_set
            all_read_set = all_read_set.union(read_name_set)
        for read_name in all_read_set:
            ref_start = ref.start
            ref_end = ref.end
            read_seq = read_name_info_dict[read_name].seq
            read_seq = read_name_info_dict[read_name].update_seq(reference_sequence=ref.reference_sequence, reference_start=ref.reference_start) if read_seq == "" else read_seq
            if not len(read_seq):
                continue
            ref_offset, alt_offset, pre_end = 0, 0, 0
            for start, end, seq in read_seq:
                if end < ref_start or start > ref_end:
                    continue
                start_offset = start - ref_start
                end_offset = end - ref_start
                if start_offset >= pre_end:
                    all_read_seq_dict[read_name] += ref.seq[pre_end:start_offset] + seq
                    pre_end = end_offset
            if pre_end < len(ref.seq):
                all_read_seq_dict[read_name] += ref.seq[pre_end:]
        read_seqs_counter = Counter(all_read_seq_dict.values())

    def extend_genotype(variants_and_genotypes, next_pos_list):

        """
        We give two iterator of two haplotype start position and update them separately, which allow haplotype extension
        more flexible when one of the haplotype has insertion or deletion, if the start position reach our provided region,
        then the extension will stop.
        """
        if next_pos_list is None or None in next_pos_list:
            pass
        if not variants_and_genotypes:
            hap1_last_pos, hap2_last_pos = next_pos_list

            rest_seq_1 = ref.query(hap1_last_pos,
                                             ref.end) if hap1_last_pos != ref.end and hap1_last_pos else ""
            rest_seq_2 = ref.query(hap2_last_pos,
                                             ref.end) if hap2_last_pos != ref.end and hap2_last_pos else ""
            yield [rest_seq_1, rest_seq_2]  # add last padding ref base
        else:
            current_variant, remaining_variants = [variants_and_genotypes[0]], variants_and_genotypes[1:]
            prefix_seqs_list, next_pos_list = find_seqs(
                current_variant, next_pos_list, ref, output_seq)
            prefix_seqs_list = list(prefix_seqs_list)

            if output_seq:
                rest_seq_1, rest_seq_2 = prefix_seqs_list
                # print(rest_seq_1, rest_seq_2)
                if len(rest_seq_1) > len(rest_seq_2):
                    rest_seq_2 += '-' * (len(rest_seq_1) - len(rest_seq_2))
                elif len(rest_seq_1) < len(rest_seq_2):
                    rest_seq_1 += '-' * (len(rest_seq_2) - len(rest_seq_1))

                prefix_seqs_list = [rest_seq_1, rest_seq_2]

            if not prefix_seqs_list or next_pos_list is None:
                pass

            for seqs in extend_genotype(remaining_variants, next_pos_list):
                yield [prefix_seqs_list[0] + seqs[0], prefix_seqs_list[1] + seqs[1]]

    def extend_genotypes(variants_and_genotypes, next_pos):

        try:
            for r in extend_genotype(variants_and_genotypes, next_pos):
                yield r
        except:
            pass

    GT = collections.namedtuple('GT', ['variant', 'genotypes'])
    if output_seq:
        variants_and_genotypes = [GT(v, g) for v, g in zip(variants, candidate_genotypes)]
        for seqs in extend_genotypes(variants_and_genotypes, [ref.start, ref.start]):
            return seqs

    genotypes_combinations = genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths, alt_dict=alt_dict, no_match=no_match_found)
    genotypes_seqs_dict = defaultdict(list)
    genotypes_list = unique_genotypes_selection(genotypes_combinations)

    for genotypes in genotypes_list:
        variants_and_genotypes = [GT(v, g) for v, g in zip(variants, genotypes)]
        for seqs in extend_genotypes(variants_and_genotypes, [ref.start, ref.start]):

            if read_seqs_counter and (seqs[0] not in read_seqs_counter or seqs[1] not in read_seqs_counter):
                continue
            genotypes_seqs_dict[frozenset(seqs)].append(genotypes)

    return genotypes_seqs_dict, read_seqs_counter

def remove_common_suffix(ref_base, alt_base):
    """
    For each haploid match, we simplify the reference base and alternative base and remove their common suffix characters.
    """

    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))  # keep at least one base
    prefix = ref_base[::-1]
    for string in alt_base:
        string = string[::-1]
        while string[:len(prefix)] != prefix and prefix:
            prefix = prefix[:len(prefix) - 1]
        if not prefix:
            break
    res_length = len(prefix)
    if res_length > min_length:
        return ref_base, alt_base
    return ref_base[:len(ref_base) - res_length], [item[:len(item) - res_length] for item in alt_base]

    return ref_base[-min_length], [item[-min_length] for item in alt_base]


def has_multi_in_truths(truths, start=None):
    for t in truths:
            if len(t.alternate_bases) > 1:
                return True
    return False

def count_combination(genotypes_combinations):
    """
    Calculate the Cartesian product required for all genotype combinations
    """
    product = 1
    for gc in genotypes_combinations:
        product *= len(gc)
    return product

def genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths=None, alt_dict=None, no_match=False, simplfy_combination=False):

    """
    Calculate genotype combination for haplotype set generation. For a locked confident variant or candidate site, we directly
    extend with its phased genotype, while for other candidates, we need to assign with a missing flag to skip the candidate, or a
    pass flag keep the genotype for further extension.
    """

    if no_match:
        [{(0, 0)}] * len(variants)
    if variant_type == 'truth':
      output = []
      for v in variants:
        if v:
          gt = tuple(v.genotype)

          is_confident_pos = (variant_dict[v.start].confident_variant and variant_dict[v.start].phased_genotype) if v.start in variant_dict else False
          output.append({tuple(variant_dict[v.start].phased_genotype)} if is_confident_pos else {(0, 0), tuple(gt), tuple(list(gt)[::-1])})
        else:
          output.append({(-1, -1)})
      return output
    elif variant_type == 'candidate':
        output = []
        has_multi_in_truth = has_multi_in_truths(truths)

        for v in variants:
            is_confident_pos = v.start in alt_dict and alt_dict[v.start].phased_genotype
            pos_in_truths = is_confident_pos and v.start in variant_dict and variant_dict[v.start].phased_genotype
            if pos_in_truths:
                output.append({tuple(alt_dict[v.start].phased_genotype)})
            elif is_confident_pos:
                if simplfy_combination:
                    output.append({tuple(alt_dict[v.start].phased_genotype)})
                else:
                    output.append({(0, 0), tuple(alt_dict[v.start].phased_genotype)})
            else:
                gt_set = set()
                for idx_1 in range(len(v.alternate_bases) + 1):
                    for idx_2 in range(len(v.alternate_bases) + 1):
                        if simplfy_combination and not has_multi_in_truth:
                                if idx_1 != 0 and idx_2 != 0 and idx_1 != idx_2:
                                    continue
                        gt_set.add((idx_1, idx_2))
                output.append(gt_set)
        # in extra mass combination, need to simplfy low possible cases:
        if  count_combination(output) > max_calculate_count:
            if simplfy_combination:
                # skip
                return [{(0, 0)}] * len(variants)
            else:
                return genotypes_combination(variants, variant_type, variant_dict, max_calculate_count, truths, alt_dict, no_match, simplfy_combination=True)
        return output


def find_seqs(variants_and_genotypes, last_pos_list, ref, output_seq=False):
    seqs = ["", ""]
    genotypes = [vg.genotypes for vg in variants_and_genotypes]
    variants = [vg.variant for vg in variants_and_genotypes]
    all_genotypes = [tuple([item[i] for item in genotypes]) for i in [0,1]]
    next_last_pos_list = [0,0]
    for idx, phased_genotype in enumerate(all_genotypes):
        current_seq, hap_end = build_seq(variants, phased_genotype, ref, last_pos_list[idx], None, output_seq=output_seq) # if overlap, merge multiple variants together
        next_last_pos_list[idx] = hap_end
        if current_seq:
            seqs[idx] = current_seq
    return seqs, next_last_pos_list


def build_seq(variants, phased_genotype, ref, pre_start, ref_end=None, output_seq=False):

    """
    Build or extend the haplotype according to provided genotype. We marked the start position iterator of each haplotype and
    update with variant alternative base.
    """

    seqs = ""
    position = pre_start
    for variant, phased in zip(variants, phased_genotype):
        if variant.start < pre_start:
            if variant.start == pre_start - 1 and phased != 0:  # this only happen when pre pos is deletion and current pos is insertion
                ref_base = variant.reference_bases
                alt_base = variant.alternate_bases[phased - 1]
                if len(alt_base) > len(ref_base): # is an insertion
                    # print ('has insertion and deletion overlap'.format(variant.start))
                    return alt_base[1:], position
            if phased != 0: # impossible # sometimes happen in true vcf
                return None, None
            else:
                return "", pre_start # do not do anything if 0 allele
        else:
            seqs += ref.query(pre_start, variant.start)

        allele = variant.reference_bases if phased == 0 else variant.alternate_bases[phased - 1]
        if phased == 0:
            allele = allele[0]
            position = variant.start + 1
            seqs += allele # only add one ref base
        else:
            ref_base = variant.reference_bases
            alt_base = variant.alternate_bases[phased-1]
            ref_base, alt_base = remove_common_suffix(ref_base, [alt_base])
            end = variant.start + len(ref_base)
            position = end
            seqs += alt_base[0].lower() if output_seq else alt_base[0]

    return seqs, position

class ReadMatch(object):
    def __init__(self, sample_ctg_info, candidates, candidate_genotypes, truths,match_seq,truth_genotypes,match_reads_count=0):

        self.sample_ctg_info = sample_ctg_info
        self.candidates = candidates
        self.truths = truths
        self.candidate_genotypes = candidate_genotypes
        self.truth_genotypes = truth_genotypes
        self.match_reads_count = match_reads_count
        self.match_seq = sorted(match_seq)
        self.truths_pos_list = [t.start for t in truths]
        self.candidates_pos_list = [c.start for c in candidates]
        self.raw_genotypes = [tuple(v.genotype) for v in self.truths]
        self.non_variants = [int(sum(cg) == 0) for cg in self.candidate_genotypes]
        self.miss_variants_count = sum([1 if sum(gt) < sum(raw_gt) else 0 for raw_gt, gt in zip(self.raw_genotypes,self.truth_genotypes)])
        self.match_variants_count = sum([1 if item < 1 else 0 for item in self.non_variants])
        self.non_variants_count = sum(self.non_variants)
        self.match_order = (self.match_reads_count, self.miss_variants_count, self.non_variants_count, self.match_variants_count)

    def match_info(self):
        can_info, truth_info = "", ""
        for can, gt in zip(self.candidates, self.candidate_genotypes):
            gt_str = '_' + '/'.join(map(str, gt)) + ' '
            can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

        for truth, gt in zip(self.truths, self.truth_genotypes):
            gt_str = '_' + '/'.join(map(str, gt)) + ' '
            truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(truth.alternate_bases) + gt_str

        extro_info = ""
        if self.match_reads_count >= -6 and self.miss_variants_count > 0 and self.match_variants_count > 0:
            extro_info = '\nthis match has few read support'
        return ('ctg_info={}, read_support,miss_variants,non_variants,match_variants={}, candidate={}, truth={} {}').format(self.sample_ctg_info,
            self.match_order,
            can_info, truth_info, extro_info)




class Read(object):
    def __init__(self, hap=0):
        self.hap = hap
        self.pos_alt_dict = defaultdict(str)
        self.start = None
        self.end = None
        self.seq = []

    def update_seq(self, reference_sequence, reference_start):
        for pos, alt_base in self.pos_alt_dict.items():
            self.start = min(self.start, pos) if self.start is not None else pos
            if alt_base[0] == 'X':
                self.seq.append((pos, pos+1, alt_base[1]))
            elif alt_base[0] == 'I':
                self.seq.append((pos, pos+1, alt_base[1:]))
            elif alt_base[0] == 'D':
                del_length = len(alt_base[1:])
                self.seq.append((pos, pos+del_length+1, reference_sequence[pos - reference_start]))
            else: #"R"
              continue
        self.end = max([item[1] for item in self.seq]) if len(self.seq) else None

        return self.seq

def decode_alt_info(cigar_count, ref_base, depth, minimum_allele_gap):
    """
    Decode the input read-level alternative information
    cigar_count: each alternative base including snp, insertion and deletion of each position
    pileup_bases: pileup bases list of each read in specific candidate position from samtools mpileup 1.10
    reference_sequence: the whole reference sequence index by contig:start-end. 0-based.
    ref_base: upper reference base for cigar calculation.
    depth: depth of candidate position for calculation.
    minimum_allele_gap: default minimum allele frequency for candidate to consider as a potential true variant for unification.
    """
    alt_type_list = []  # SNP I D
    seqs = cigar_count.split(' ')
    seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}
    max_del_cigar = ""
    del_list = []
    ref_represatation = ref_base
    alt_list = sorted(list(seq_alt_bases_dict.items()), key=lambda x: x[1], reverse=True)

    # select maximum 2 variants type
    seq_insertion_bases_list = alt_list[:2]
    for alt_type, count in seq_insertion_bases_list:
        count = int(count)
        if count / float(depth) < minimum_allele_gap:
            continue
        if alt_type[0] == 'X':
            alt_type_list.append(alt_type[1])
        elif alt_type[0] == 'I':
            alt_type_list.append(alt_type[1:])
        elif alt_type[0] == 'D':
            if len(alt_type[1:]) > len(max_del_cigar):
                max_del_cigar = alt_type[1:]
            del_list.append(ref_base + alt_type[1:])
    new_del_list = []
    if len(max_del_cigar):
        ref_represatation = ref_base + max_del_cigar
        alt_type_list = [item + max_del_cigar for item in alt_type_list]
        for idx, item in enumerate(del_list):
            start_pos = len(item[1:])
            append_del_bases = max_del_cigar[start_pos:]
            new_del_list.append(
                ref_base + append_del_bases)  # ACG-> A, ACGTT -> A, max_del_cigar is CGTT, represent ACG-> A to ACGTT->ATT
    alt_base_list = alt_type_list + new_del_list
    return ref_represatation, alt_base_list, alt_list

def has_variant_suport(ref_base, alt_base, pos, alt_dict):
    """
    ref_base: reference base of the true varaint.
    alt_base: alternative base of the true variant
    pos: pos: candidate position for unification.
    alt_dict: dictionary (pos: pos info) which keep position level candidate reference base and alternative base information.
    return the alternative index of each candidate site if found match
    """

    alt_index = -1
    if pos not in alt_dict or not len(alt_dict[pos]):
        return alt_index
    cigar_count = alt_dict[pos]
    if len(ref_base) == 1 and len(alt_base) == 1:  # Snp
        if alt_base in cigar_count:
            alt_index = cigar_count[0][alt_base]
    elif len(ref_base) > len(alt_base):  # D
        del_base = ref_base[1:len(ref_base) - len(alt_base) + 1]
        if del_base in cigar_count:
            alt_index = cigar_count[1][alt_base]
    elif len(ref_base) < len(alt_base):  # I
        ins_base = alt_base[1:len(alt_base) - len(ref_base) + 1]
        if ins_base in cigar_count:
            alt_index = cigar_count[2][alt_base]
    return alt_index


def get_ref(ref_fn, ctg_name):
    refernce_sequences = []
    samtools_faidx_process = subprocess_popen(
        shlex.split("samtools faidx {} {}".format(ref_fn, ctg_name))
    )
    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    reference_sequence = "".join(refernce_sequences[1:])

    reference_sequence = reference_sequence.upper()
    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    reference_start = 1
    return reference_sequence, reference_start


def get_genotype(genotype):
    g1, g2 = genotype
    min_gt = min(int(g1), int(g2))
    max_gt = max(int(g1), int(g2))
    return str(min_gt) + '/' + str(max_gt)

def lock_variant(variant, truth):
    """
    Find potential locked true variant and candidate site if we consider it as a confident site, there are only exactly one
    candidate match with the true variant in a specific site, and we can further take the candidate into consideration with
    a match flag.
    """
    if truth is None:
        return None, False
    variant_ref_base = variant.reference_bases
    variant_alt_base = variant.alternate_bases
    truth_ref_base = truth.reference_bases
    truth_alt_base = truth.alternate_bases

    if len(variant_alt_base) != len(truth_alt_base):
        return None, False

    tmp_alt_list = []
    for ab in variant_alt_base:
        ref_base1, alt_base1 = remove_common_suffix(variant_ref_base, [ab])
        tmp_alt_list.append((alt_base1[0], ref_base1))
    match_index = [-1] * len(truth_alt_base)

    for t_idx, ab in enumerate(truth_alt_base):
        ref_base1, alt_base1 = remove_common_suffix(truth_ref_base, [ab])
        for a_idx, (alt_base, ref_base) in enumerate(tmp_alt_list):
            if alt_base1[0] == alt_base and ref_base1 == ref_base:
                match_index[t_idx] = a_idx
    match = sum([1 if item >= 0 else 0 for item in match_index]) == len(truth_alt_base) # can find all alt_base
    return match_index, match

def decode_variant(variant, reference_base):
    if variant == 'R':
        return 'R', 'R'
    if variant[0] == 'X':
        return reference_base, variant[1]
    elif variant[0] == 'I':
        return reference_base,  variant[1:]
    elif variant[0] == 'D':
        return reference_base + variant[1:], reference_base

def update_variant_hap_dict(alt_dict, pos, reference_sequence, reference_start, is_variant_confident, variant_dict, allele_gap, platform):
    """
    For a phased alignment, the candidates are easier to lock as confident if the signal exists strongly in one side and have confident
    match with true variant.
    """
    phased_genotype = [-1,-1]
    reference_base = reference_sequence[pos - reference_start]
    variant_hap_dict = alt_dict[pos].variant_hap_dict
    if not len(variant_hap_dict):
        return None
    hap_count_dict = alt_dict[pos].hap_count_dict
    variant_ref_base = alt_dict[pos].reference_bases
    variant_alt_base = alt_dict[pos].alternate_bases
    for variant, hap_dict in variant_hap_dict.items():
        if variant in '*#':
            continue
        hap_0 = hap_dict[0] if 0 in hap_dict else 0
        # for illumina unification, phased information contributes less, we safely denote with reference allele gap
        if platform == 'ilmn':
            hap_total_af = (hap_0) / float(sum(list(hap_count_dict.values())))
            if variant not in 'R*' and hap_total_af > reference_allele_gap and is_variant_confident and variant_dict[pos].genotype == [1, 1]:
                phased_genotype = [1, 1]
                return phased_genotype
            if -1 in phased_genotype:
                return None
            return phased_genotype

        hap_1 = hap_dict[1] if 1 in hap_dict else 0
        hap_2 = hap_dict[2] if 2 in hap_dict else 0
        hap_0 = hap_0 if hap_0 > 3 else 0
        hap_1 = hap_1 if hap_1 > 3 else 0
        hap_2 = hap_2 if hap_2 > 3 else 0

        hap_total_af = (hap_0 + hap_1 + hap_2) / float(sum(list(hap_count_dict.values())))
        if variant not in 'R*#' and hap_total_af > 1 - allele_gap / 2 and is_variant_confident and variant_dict[pos].genotype== [1, 1]:
            phased_genotype = [1, 1]
            return phased_genotype
        hap_1_af = hap_1 / float(hap_count_dict[1]) if 1 in hap_count_dict else 0.0
        hap_2_af = hap_2 / float(hap_count_dict[2]) if 2 in hap_count_dict else 0.0

        if variant == 'R':
            if hap_1_af >= 1 - allele_gap:
                phased_genotype[0] = 0
            if hap_2_af >= 1 - allele_gap:
                phased_genotype[1] = 0
            continue
        ref_base, alt_base = decode_variant(variant, reference_base)

        alt_index = -1
        for ab_idx, ab in enumerate(variant_alt_base):
            ref_base1, alt_base1 = remove_common_suffix(variant_ref_base, [ab])
            if alt_base1[0] == alt_base and ref_base1 == ref_base:
                alt_index = ab_idx
                break
        if alt_index == -1:
            continue
        if hap_1_af >= 1 - allele_gap * 2:
            phased_genotype[0] = alt_index + 1
        if hap_2_af >= 1 - allele_gap * 2:
            phased_genotype[1] = alt_index + 1
    del alt_dict[pos].variant_hap_dict
    del alt_dict[pos].hap_count_dict
    if -1 in phased_genotype:
        return None
    return phased_genotype


def match_alt_base(alt_list, ref_base, alt_base):
    if not len(alt_list) or (len(alt_list) == 1 and 'R' in alt_list):
        return False
    alt_set = set([item[0] for item in alt_list])

    for ab in alt_base:
        rb, ab = remove_common_suffix(ref_base, [ab])
        if len(rb) == len(ab[0]): #snp
            ab = 'X' + ab[0]
            if ab in alt_set:
                return True
        elif len(rb) < len(ab[0]): # insertion
            ab = 'I' + ab[0]
            if ab in alt_set:
                return True
        elif len(rb) > len(ab[0]):
            ab = 'D' + rb[1:]
            if ab in alt_set:
                return True
    return False

def check_confident_match(candidates, truths):

    """
    Double check whether the candidate site match the representation in the truth variant site in reference base,
    alternative base, genotype and start position.
    """

    if len(candidates) != len(truths):
        return False
    all_candidate_positions = set([c.start for c in candidates])
    for truth in truths:
        if truth.start not in all_candidate_positions:
            return False
        for candidate in candidates:
            if candidate.start == truth.start:
                if candidate.reference_bases != truth.reference_bases or candidate.alternate_bases != truth.alternate_bases:
                    return False
    return True

INFO = collections.namedtuple('INFO', ['start', 'type', 'variant'])
class SplitVariants(object):
    def __init__(self, truths,
                 max_candidates_distance,
                 max_calculate_count,
                 partition_size,
                 variant_dict=None):

        self.max_candidates_distance = max_candidates_distance
        self.max_calculate_count = max_calculate_count
        self.partition_size = partition_size
        self.truths = truths
        self.product_count = 1
        self.variant_dict = variant_dict
        self.truths_pos_list = list(sorted([v.start for v in truths]))
        self.truths_pos_set = set(self.truths_pos_list)
        self.truth_index = 0
        self.truth_start = self.truths_pos_list[0]
        self.partition = []

    def all_genotypes_combination(self, variant):
        """
        Enumerate true variant site and candidate site genotype combination and find read. For a phased confident site, we
        only enumerate one genotype with confident flag.
        """

        if variant.type == 'candidate':
            num_ref_and_alts = len(variant.variant.alternate_bases)
            candidate = variant.variant
            if candidate.phased_genotype:
                return 1
            elif candidate.confident_variant:
                return 2
            return (num_ref_and_alts + 1) * num_ref_and_alts / 2
        else:
            if variant.variant.phased_genotype is not None:
                return 1
            return len(variant.variant.alternate_bases)

    def match_max_candidate_distance(self, variants, new_count):
        if not len(self.partition):
            return True
        n_of_type = sum(1 for g in self.partition if g.type == variants.type)
        if new_count >= max_calculate_count or n_of_type >= self.partition_size:
            if new_count >= max_calculate_count:
                print('{} exceed max calculation count'.format(new_count))
            return False
        else:
            for g in self.partition:
                if variants.variant.start - g.variant.end + 1 > self.max_candidates_distance:
                    return False

            last_par = self.partition[-1].variant.end
            if variants.variant.start - last_par + 1 > extend_bp:
                return False
            return True

    def merge_truth(self, candidate):
        all_variants = []
        candidate_start = candidate.start
        if self.truth_start is None or candidate_start < self.truth_start:
            all_variants.append(INFO(candidate.start, 'candidate', candidate))
        elif self.truth_start is not None:
            all_variants.append(INFO(candidate.start, 'truth', self.truths[self.truth_index]))
            while candidate_start >= self.truth_start:
                truth = self.truths[self.truth_index]
                all_variants.append(INFO(truth.start, 'truth', truth))
                if self.truth_index == len(self.truths_pos_list) - 1:
                    all_variants.append(INFO(candidate.start, 'candidate', candidate))
                    break
                else:
                    self.truth_index += 1
                    self.truth_start = self.truths_pos_list[self.truth_index]
        return all_variants

    def split_truths_and_candidates(self, partition):
        candidate_partitions = []
        truth_partitions = []
        for p in partition:
            if p.type == 'candidate':
                candidate_partitions.append(p.variant)
            elif p.type == 'truth':
                truth_partitions.append(p.variant)
        return [candidate_partitions, truth_partitions]

    def add_variant(self, truths_reader, candidates_reader):

        for variants in heapq.merge(truths_reader, candidates_reader):
            new_count = self.product_count * self.all_genotypes_combination(
                variants)

            if self.match_max_candidate_distance(variants,
                                            new_count):
                self.partition.append(variants)
                self.product_count = new_count
            else:
                if variants.start == self.partition[-1].start and variants.type != self.partition[
                    -1].type:  #
                    # add same truths or variants together and add at least one nearby candidate
                    self.partition.append(variants)
                    # if sv_idx < len(sorted_variants) - 1 and sorted_variants[sv_idx + 1].start not in truths_pos_set and \
                    #         sorted_variants[sv_idx + 1].start - variants.start <= extend_bp:
                    #     partition.append(sorted_variants[sv_idx + 1])
                    partition = self.partition
                    self.partition = []
                    self.product_count = 1
                    yield self.split_truths_and_candidates(partition)
                else:
                    partition = self.partition
                    self.partition = [variants]
                    self.product_count = self.all_genotypes_combination(variants)
                    yield self.split_truths_and_candidates(partition)
        if len(self.partition):
            yield self.split_truths_and_candidates(self.partition)

class RepresentationUnification(object):

    def __init__(self,
                 sample_name,
                 contig_name,
                 reference_sequence,
                 reference_start,
                 partition_size,
                 max_candidates_distance,
                 max_calculate_count,
                 subsample_ratio):

        self.sample_name = sample_name
        self.contig_name = contig_name
        self.subsample_ratio = subsample_ratio
        self.sample_ctg_info = '_'.join([sample_name, str(subsample_ratio), contig_name])
        self.partition_size = partition_size
        self.max_candidates_distance = max_candidates_distance
        self.max_calculate_count = max_calculate_count
        self.reference_sequence = reference_sequence
        self.reference_start = reference_start


    def get_reference_seq(self, candidates, true_variants, bufsize=50):
        all_variants = candidates + true_variants
        start = min(x.start for x in all_variants)
        end = max(x.end for x in all_variants)

        ref_bases = self.reference_sequence[start - self.reference_start - 1:end + bufsize - self.reference_start]
        return Reference(seq=ref_bases,
                         start=start - 1,
                         reference_sequence=self.reference_sequence,
                         reference_start=self.reference_start)

    def find_match_pairs(self, candidates, truths, ref, variant_dict, read_name_info_dict=None, alt_dict=None):
        no_match_found = len(candidates) == 0 or len(truths) == 0

        if test_pos:
            print('variants:')
            for v in candidates:
                print (v.start, v.reference_bases, v.alternate_bases, v.phased_genotype)
            print ('truths:')
            for v in truths:
                print (v.start, v.reference_bases, v.alternate_bases, v.phased_genotype, v.confident_variant)
            print('\n')

        if no_match_found:
            can_info, truth_info = "", ""
            for can in candidates:
                gt = can.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

            for truth in truths:
                gt = truth.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(
                    truth.alternate_bases) + gt_str

            print ('[INFO] Missing match: ctg_info={}, read_support,miss_variants,non_variants,match_variants=(0, {}, {}, 0), candidate={}, truth={}'.format(self.sample_ctg_info,len(truths), len(candidates),can_info, truth_info))
            return None

        confident_match = check_confident_match(candidates, truths)
        if confident_match:
            can_info, truth_info = "", ""
            for can in candidates:
                gt = can.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                can_info += str(can.start) + '-' + can.reference_bases + '->' + '-'.join(can.alternate_bases) + gt_str

            for truth in truths:
                gt = truth.genotype
                gt_str = '_' + '/'.join(map(str, gt)) + ' '
                truth_info += str(truth.start) + '-' + truth.reference_bases + '->' + '-'.join(
                    truth.alternate_bases) + gt_str

            print (
                '[INFO] Found confident match: ctg_info={}, read_support,miss_variants,non_variants,match_variants=(None, 0, 0, {}), candidate={}, truth={}'.format(
                    self.sample_ctg_info, len(truths), can_info, truth_info))

            match_genotype = [tuple(v.genotype) for v in truths]

            return ReadMatch(
                sample_ctg_info=self.sample_ctg_info,
                candidates=truths,
                candidate_genotypes=match_genotype,
                truths=truths,
                truth_genotypes=match_genotype,
                match_seq=ref.seq,
                match_reads_count=100)

        truths_candidate_gentoypes = genotypes_combination(truths, 'truth', variant_dict, max_calculate_count,
                                                       truths, alt_dict=alt_dict, no_match=no_match_found)
        candidates_candidate_gentoypes = genotypes_combination(candidates, 'candidate', variant_dict,
                                                                              max_calculate_count,
                                                                              truths, alt_dict=alt_dict,
                                                                              no_match=no_match_found)

        truths_genotypes_list = unique_genotypes_selection(truths_candidate_gentoypes)
        candidates_genotypes_list = unique_genotypes_selection(candidates_candidate_gentoypes)

        # print (len(truths_genotypes_list) * len(candidates_genotypes_list))
        if len(truths_genotypes_list) * len(candidates_genotypes_list) > self.max_calculate_count:
            return None
        variant_seqs, read_seqs_counter = find_read_support(
            variants=candidates,
            ref=ref,
            variant_type='candidate',
            max_calculate_count=self.max_calculate_count,
            variant_dict=variant_dict,
            truths=truths,
            read_name_info_dict=read_name_info_dict,
            alt_dict=alt_dict,
            no_match_found=no_match_found)

        truth_seqs, _ = find_read_support(
            variants=truths,
            ref=ref,
            variant_type='truth',
            max_calculate_count=self.max_calculate_count,
            variant_dict=variant_dict,
            truths=truths,
            read_seqs_counter=read_seqs_counter,
            read_name_info_dict=read_name_info_dict,
            alt_dict=alt_dict,
            no_match_found=no_match_found)

        matches = []
        for variant_seq, variant_genotypes in variant_seqs.items():
            if variant_seq not in truth_seqs:
                continue
            truth_seq = truth_seqs[variant_seq]
            for variant_genotype in variant_genotypes:
                match_reads_count = -sum([read_seqs_counter[seq] if read_seqs_counter and seq in read_seqs_counter else 0 for seq in variant_seq ]) # more match reads and negative count is better if smaller
                matches.append(ReadMatch(
                        sample_ctg_info=self.sample_ctg_info,
                        candidates=candidates,
                        candidate_genotypes=variant_genotype,
                        truths=truths,
                        truth_genotypes=truth_seq[0],
                        match_seq=variant_seq,
                        match_reads_count=match_reads_count))
        if not matches:
            return None
        else:
            best_matches = sorted(matches, key=lambda x: x.match_order)[0]
            # print ('[INFO] Found match case:', best_matches.match_info())
            output_seq = False
            if output_seq:
                candidate_genotypes = best_matches.candidate_genotypes
                truth_genotypes = best_matches.truth_genotypes
                candidate_seqs = find_read_support(candidates,ref, 'candidate', candidate_genotypes=candidate_genotypes, max_calculate_count=self.max_calculate_count,variant_dict=variant_dict, truths=truths, read_name_info_dict=read_name_info_dict,  alt_dict=alt_dict, output_seq=True)
                truths_seqs = find_read_support(truths,ref, 'truth', candidate_genotypes=truth_genotypes, max_calculate_count=self.max_calculate_count,variant_dict=variant_dict, truths=truths, read_name_info_dict=read_name_info_dict,  alt_dict=alt_dict, output_seq=True)
                # print (set(candidate_seqs) == set(truths_seqs))
                # update_candidate_seqs = set([item.replace('-', '').upper() for item in candidate_seqs])
                # update_truths_seqs = set([item.replace('-', '').upper() for item in truths_seqs])
                # print('[INFO] match? ', set(update_candidate_seqs) == set(update_truths_seqs))
                sorted_truths_seqs = truths_seqs if truths_seqs[0].replace('-', '').upper() ==  candidate_seqs[0].replace('-', '').upper() else truths_seqs[::-1]
                output_info = [self.sample_name, self.contig_name, str(self.subsample_ratio), str(ref.start), str(ref.end), ref.seq, ','.join(candidate_seqs), ','.join(sorted_truths_seqs) ]
                print ('[INFO] Found match case:', best_matches.match_info() + '~' + '~'.join(output_info) + '~'+ ','.join([str(t.pos) for t in candidates]))
            else:
                print ('[INFO] Found match case:', best_matches.match_info())

            return best_matches


    def unify_label(self, variants, truths, ctg_start, ctg_end, all_pos, variant_dict,
                  rescue_dict=None, output_vcf_fn=None, read_name_info_dict=None, alt_dict=None):


        all_candidates, all_truths = variants, truths
        ref = self.get_reference_seq(all_candidates, all_truths)
        match_pairs = self.find_match_pairs(candidates=all_candidates,
                                            truths=all_truths,
                                            ref=ref,
                                            variant_dict=variant_dict,
                                            read_name_info_dict=read_name_info_dict,
                                            alt_dict=alt_dict)

        if match_pairs is None:
            if not len(truths):
                return
            # double check to rescue true variants
            for truth in all_truths:
                pos = truth.start
            # add missing low-confident tp position
                if not (pos >= ctg_start and pos < ctg_end):
                    continue
                if pos in alt_dict and pos in variant_dict:
                    ref_base = variant_dict[pos].reference_bases
                    alt_base = variant_dict[pos].alternate_bases
                    alt_list = alt_dict[pos].alt_list
                    if not match_alt_base(alt_list, ref_base, alt_base):
                        print('[INFO] {} {} miss and has no cigar support'.format(self.sample_ctg_info, pos))
                        continue
                    print('[INFO] {} {} miss by match, append to vcf'.format(self.sample_ctg_info, pos))
                    if pos in all_pos or pos in rescue_dict:
                        continue
                    ref_base = variant_dict[pos].reference_bases
                    variant = ','.join(variant_dict[pos].alternate_bases)
                    genotype_string = '/'.join(map(str, variant_dict[pos].genotype))
                    # For efficiency, we currently only compute reference base, altnertive base and genotype from GetTruth.py
                    rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                            self.contig_name,
                            pos,
                            ref_base,
                            variant,
                            10,
                            'PASS',
                            'R',
                            genotype_string,
                            10,
                            10,
                            0.5)
                else:
                    print('[INFO] {} {} miss and no variant support'.format(self.sample_ctg_info, pos))
            return
        for idx, (candidate, candidate_genotypes) in enumerate(
                zip(match_pairs.candidates, match_pairs.candidate_genotypes)):
            pos = candidate.start

            have_miss_variants = True if sum([1 for gt in match_pairs.truth_genotypes if sum(gt) == 0]) else False

            # append a position into rescue queue if it was missed by the unification
            if sum(candidate_genotypes) == 0 and pos not in rescue_dict and have_miss_variants and pos not in variant_dict and pos in alt_dict and alt_dict[pos].phased_genotype:
                genotype_string = '/'.join(map(str, alt_dict[pos].phased_genotype))
                variant = ','.join(candidate.alternate_bases)
                ref_base = candidate.reference_bases
                rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                    self.contig_name, pos, ref_base, variant, 10, 'PASS', 'R', genotype_string, 10, 10, 0.5)
                continue
            if sum(candidate_genotypes) == 0:
                continue
            if not len(candidate.alternate_bases):
                continue
            g1, g2 = candidate_genotypes
            variant = set()
            ref_base = candidate.reference_bases
            if g1 != 0:
                variant.add(candidate.alternate_bases[g1 - 1])
            if g2 != 0:
                variant.add(candidate.alternate_bases[g2 - 1])
            if g1 == 0 or g2 == 0:
                genotype_string = '0/1'
            elif g1 == g2:
                genotype_string = '1/1'
            elif g1 != g2:
                genotype_string = '1/2'
            ref_base, variant = remove_common_suffix(ref_base, list(variant))
            variant = ','.join(variant)
            if candidate.start in all_pos:
                continue
            all_pos.add(pos)

            if output_vcf_fn is not None:
                # For efficiency, we only compute reference base, altnertive base and genotype for GetTruth.py currently
                print("%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                    self.contig_name, candidate.start, ref_base, variant, 10, 'PASS', 'U', genotype_string, 10, 10, 0.5), file=output_vcf_fn)
                if pos in rescue_dict:
                    del rescue_dict[pos]
        for idx, (pos, raw_genotype, truth_genotype) in enumerate(
                zip(match_pairs.truths_pos_list, match_pairs.raw_genotypes,
                    match_pairs.truth_genotypes)):
            if not (pos >= ctg_start and pos < ctg_end):
                continue

            if truth_genotype == (0, 0) and sum(raw_genotype) > 0:# miss genoytpe
                if pos in alt_dict and pos in variant_dict:
                    ref_base = variant_dict[pos].reference_bases
                    alt_base = variant_dict[pos].alternate_bases
                    alt_list = alt_dict[pos].alt_list
                    if not match_alt_base(alt_list, ref_base, alt_base):
                        print('{} {} miss and has no cigar support'.format(self.sample_ctg_info, pos))
                        continue
                    print('{} {} miss by match, append to vcf'.format(self.sample_ctg_info, pos))
                    if pos in all_pos:
                        continue
                    all_pos.add(pos)

                    ref_base = variant_dict[pos].reference_bases
                    variant = ','.join(variant_dict[pos].alternate_bases)
                    genotype_string = '/'.join(map(str, variant_dict[pos].genotype))

                    if output_vcf_fn is not None:
                        rescue_dict[pos] = "%s\t%d\t.\t%s\t%s\t%d\t%s\t%s\tGT:GQ:DP:AF\t%s:%d:%d:%.4f" % (
                            self.contig_name, pos, ref_base, variant, 10, 'PASS', 'R', genotype_string, 10, 10, 0.5)


def find_chunk_id(center_pos, chunk_num, ctg_name, fai_fn):
        contig_length = 0
        with open(fai_fn, 'r') as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")

                contig_name = columns[0]
                if contig_name != ctg_name:
                    continue
                contig_length = int(columns[1])
        chunk_size = contig_length // chunk_num + 1 if contig_length % chunk_num else contig_length // chunk_num
        for chunk_id in range(chunk_num):
            ctg_start = chunk_size * chunk_id  # 0-base to 1-base
            ctg_end = ctg_start + chunk_size
            if center_pos >= ctg_start and center_pos < ctg_end:
                return chunk_id
        return None


def UnifyRepresentation(args):

    """
    Representation Unification algorithm main function, this algorithm aims to unify variant representation
    between training material and true variant set.
    All candidate sites with sufficient read support and over a certain allele frequency were selected, of which the
    same variant information with true set was locked as confident candidate sites. Secondly, for each remaining
    candidate site and true variant site, a matching or missing flag was assigned. We build  haplotype pairs based on
    all possible flag combinations. For each fully matching pair, we use the match pair with the most support reads as
    the final best match pair. For the remaining unmatched sites, we decrease allele frequency to further seek remaining
    candidate sites. We rescue the candidate sites matching true variant type. In the end, the unified VCF consists of
    locked variant sites, unified variant sites, and rescued variant sites.
    """
    sample_name = args.sampleName
    vcf_fn = args.vcf_fn
    candidates_fn = args.candidates_fn
    contig_name = args.ctgName
    ctg_start = args.ctgStart
    ctg_end = args.ctgEnd
    bed_fn = args.bed_fn
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num
    fasta_file_path = args.ref_fn
    is_confident_bed_file_given = bed_fn is not None
    partition_size = args.partition_size
    minimum_allele_gap = args.minimum_allele_gap
    max_candidates_distance = args.max_candidates_distance
    global max_calculate_count
    max_calculate_count = args.max_calculate_count
    subsample_ratio = args.subsample_ratio
    platform = args.platform

    global test_pos
    test_pos = None
    test_pos = 29703861
    if args.test_pos and test_pos:
        platform ="hifi"
        sample_name = 'hg002'
        fasta_file_path = '/mnt/bal36/zxzheng/testData/ont/data/GRCh38_no_alt_analysis_set.fasta'
        subsample_ratio = 1000
        bed_fn = "/autofs/bal33/zxzheng/data/vcf/hg38/{}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed".format(sample_name.upper())
        is_confident_bed_file_given = True
        contig_name = 'chr6'
        # var_read_fn = '/autofs/bal33/zxzheng/pb/hap_labler/best_0/alt_info/hg002_1000_{}'.format(contig_name[3:])
        vcf_fn = '/autofs/bal33/zxzheng/data/var/{}/var_{}'.format(sample_name, contig_name[3:])
        args.output_vcf_fn = '/autofs/bal33/zxzheng/use_test/tmp_vcf'
        ctg_start = test_pos - 1000
        ctg_end = test_pos + 1000
        chunk_num = 50
        fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')
        chunk_id = find_chunk_id(center_pos=test_pos, chunk_num=50, ctg_name=contig_name, fai_fn=fai_fn)

        candidates_fn = "/mnt/bal36/zxzheng/ru/{}/{}_{}/candidates/{}_{}_{}_{}".format(platform, sample_name, contig_name[3:], sample_name, subsample_ratio, contig_name[3:], chunk_id+1)
        print ("candidate file {}".format(candidates_fn))
        chunk_id = None
    else:
        test_pos = None
    alt_dict = defaultdict()
    read_name_info_dict = defaultdict(Read)

    fai_fn = file_path_from(fasta_file_path, suffix=".fai", exit_on_not_found=True, sep='.')

    if chunk_id is not None:

        """
        Whole genome calling option, acquire contig start end position from reference fasta index(.fai), then split the
        reference according to chunk id and total chunk numbers.
        """
        contig_length = 0
        with open(fai_fn, 'r') as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")

                ctg_name = columns[0]
                if contig_name != ctg_name:
                    continue
                contig_length = int(columns[1])
        chunk_size = contig_length // chunk_num + 1 if contig_length % chunk_num else contig_length // chunk_num
        ctg_start = chunk_size * chunk_id  # 0-base to 1-base
        ctg_end = ctg_start + chunk_size

    is_ctg_name_given = contig_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None
    ref_regions = []
    reference_start = 1
    if is_ctg_range_given:
        reference_start, reference_end = ctg_start - param.expandReferenceRegion, ctg_end + param.expandReferenceRegion
        reference_start = 1 if reference_start < 1 else reference_start
        ref_regions.append(region_from(ctg_name=contig_name, ctg_start=reference_start, ctg_end=reference_end))
    elif is_ctg_name_given:
        ref_regions.append(region_from(ctg_name=contig_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command='samtools',
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )

    if is_confident_bed_file_given:
        tree = bed_tree_from(bed_fn, contig_name=contig_name)

    vcf_reader = VcfReader(vcf_fn=vcf_fn, ctg_name=contig_name, ctg_start=ctg_start-extended_window_size, ctg_end=ctg_end+extended_window_size, is_var_format=True)
    vcf_reader.read_vcf()
    variant_dict = vcf_reader.variant_dict

    # truths = sorted([(item, variant_dict[item]) for item in variant_dict.keys() if
    #                   is_region_in(
    #                      tree=tree,
    #                      contig_name=contig_name,
    #                      region_start=item-2,
    #                      region_end=variant_dict[item].end + 2)], key=lambda x: x[0])
    truths = sorted([(item, variant_dict[item]) for item in variant_dict.keys()], key=lambda x: x[0])
    truths = [item[1] for item in truths]

    if not len(truths) or not len(variant_dict):
        return

    rescue_dict = defaultdict()
    all_pos = set()

    output_vcf_fn = None
    if args.output_vcf_fn:
      output_vcf_fn = open(args.output_vcf_fn, "w")
      header = output_header(reference_file_path=fasta_file_path, sample_name=sample_name)
      output_vcf_fn.write(header)

    # For pacbio hifi platform, we select larger candidate distance for better performance
    if platform == 'hifi':
        max_candidates_distance = 200
        partition_size = 20

    RU = RepresentationUnification(
        sample_name=sample_name,
        contig_name=contig_name,
        reference_sequence=reference_sequence,
        reference_start=reference_start,
        partition_size=partition_size,
        max_candidates_distance=max_candidates_distance,
        max_calculate_count=max_calculate_count,
        subsample_ratio=subsample_ratio)

    SP = SplitVariants(truths=truths,
                 max_candidates_distance=max_candidates_distance,
                 max_calculate_count=max_calculate_count,
                 partition_size=partition_size,
                 variant_dict=variant_dict)

    def extract_candidates_generator_from():
        unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (candidates_fn)))
        for row in unzip_process.stdout:
            row = row.strip().split('\t', maxsplit=2)  # ['chr_pos', 'depth', 'cigar_count']
            chr_pos, depth, var_read_json = row[:3]
            ctg_name, pos = chr_pos.split()
            pos, depth = int(pos), int(depth)
            if contig_name and contig_name != ctg_name:
                continue
            if pos < ctg_start - extended_window_size or pos > ctg_end + extended_window_size:
                continue

            if test_pos and (pos < test_pos - 1000 or pos > test_pos + 1000):
                continue
            var_read_dict = json.loads(var_read_json)
            if not len(var_read_dict):
                continue

            cigar_count = ' '.join([' '.join([item, str(len(var_read_dict[item].split(' ')))]) for item in var_read_dict.keys()])
            ref_base = reference_sequence[pos - reference_start]
            pos_in_truths = pos in variant_dict
            ref_base, alt_base, alt_list = decode_alt_info(cigar_count=cigar_count,
                                                                   ref_base=ref_base,
                                                                   depth=depth,
                                                                   minimum_allele_gap=minimum_allele_gap)

            alt_dict[pos] = Position(pos=pos,
                                    ref_base=ref_base,
                                     alt_base=alt_base,
                                     genotype1=-1,
                                     genotype2=-1,
                                     candidate=True,
                                     depth=depth,
                                     alt_list=alt_list)

            for variant, read_str in var_read_dict.items():
                read_list = read_str.split(' ')
                for read_name in read_list:
                    read_name, hap = read_name[:-1], read_name[-1]
                    if read_name not in read_name_info_dict or read_name_info_dict[read_name].hap == 0 and hap != 0:
                        read_name_info_dict[read_name].hap = int(hap)

                    read_hap = read_name_info_dict[read_name].hap if read_name in read_name_info_dict else 0
                    read_name_info_dict[read_name].pos_alt_dict[pos] = variant
                    if read_hap in alt_dict[pos].variant_hap_dict[variant]:
                        alt_dict[pos].variant_hap_dict[variant][read_hap] += 1
                    else:
                        alt_dict[pos].variant_hap_dict[variant][read_hap] = 1
                    alt_dict[pos].hap_count_dict[read_hap] += 1
                    alt_dict[pos].read_name_set.add(read_name)

            match_index, is_variant_confident = lock_variant(alt_dict[pos], variant_dict[pos] if pos_in_truths else None)
            if is_variant_confident:
                variant_dict[pos].confident_variant = match_index
            alt_dict[pos].phased_genotype = update_variant_hap_dict(alt_dict=alt_dict,
                                                                    pos=pos,
                                                                    reference_sequence=reference_sequence,
                                                                    reference_start=reference_start,
                                                                    is_variant_confident=is_variant_confident,
                                                                    variant_dict=variant_dict,
                                                                    allele_gap=minimum_allele_gap,
                                                                    platform=platform)

            # lock the candidate if it has meet the phased_genotype requirement and have a exactly one match true variant site
            if alt_dict[pos].phased_genotype and pos_in_truths and is_variant_confident:
                if alt_dict[pos].phased_genotype.count(0) != variant_dict[pos].genotype.count(0) or (sum(variant_dict[pos].genotype) == 3 and sum(alt_dict[pos].phased_genotype) != 3):
                    # skip wrong genotype
                    alt_dict[pos].phased_genotype = None
                variant_dict[pos].reference_bases = alt_dict[pos].reference_bases
                variant_dict[pos].alternate_bases = alt_dict[pos].alternate_bases
                variant_dict[pos].phased_genotype = alt_dict[pos].phased_genotype

            yield INFO(pos, 'candidate', alt_dict[pos])

    candidates_reader = extract_candidates_generator_from()
    truths_reader = [INFO(t.start, 'truth', t) for t in truths]

    for partition in SP.add_variant(truths_reader, candidates_reader):
        # print (len(partition[0]), len(partition[1]))
        if partition is None:
            continue
        variants, truths = partition
        if test_pos and (test_pos not in set([v.start for v in variants])) and (test_pos not in set([t.start for t in truths])):
            continue
        RU.unify_label(variants=variants,
                         truths=truths,
                         ctg_start=ctg_start,
                         ctg_end=ctg_end,
                         all_pos=all_pos,
                         variant_dict=variant_dict,
                         rescue_dict=rescue_dict,
                         output_vcf_fn=output_vcf_fn,
                         read_name_info_dict=read_name_info_dict,
                         alt_dict=alt_dict)

        # empty dict
        min_start = min(min([v.start for v in variants] if len(variants) else [0]), min([t.start for t in truths] if len(truths) else [0]))
        for read_name in list(read_name_info_dict.keys()):
            read = read_name_info_dict[read_name]
            if read.end is not None and read.end < min_start:
                del read_name_info_dict[read_name]

        for pos in list(sorted(alt_dict.keys())):
            if pos < min_start:
                del alt_dict[pos]
            else:
                break

    if not len(rescue_dict):
        return
    if output_vcf_fn is not None:
        for pos, vcf_info in rescue_dict.items():
            print(vcf_info, file=output_vcf_fn)
        output_vcf_fn.close()

    if os.path.exists(args.output_vcf_fn):
        for row in open(args.output_vcf_fn, 'r'):
            if row[0] != '#':
                return
        os.remove(args.output_vcf_fn)
        print("[INFO] No vcf output for file {}, remove empty file".format(args.output_vcf_fn))

def main():
    parser = ArgumentParser(description="Representation unification for candidate site and true variant")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Truth variants or called variants")

    parser.add_argument('--ref_fn', type=str, default="ref.fa",
                        help="Reference fasta file input, default: %(default)s")

    parser.add_argument('--candidates_fn', type=str, default=None,
                        help="All read-level candidate details of provided contig, default: %(default)s")

    parser.add_argument('--output_vcf_fn', type=str, default=None,
                        help="VCF output filename or stdout if not set,default: %(default)s")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of the sequence to be processed")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    # options for advanced users
    parser.add_argument('--max_candidates_distance', type=int, default=100,
                        help="EXPERIMENTAL: Maximum distance between subsequent variants within a group")

    parser.add_argument('--max_calculate_count', type=int, default=10000,
                        help="EXPERIMENTAL: Maximum calculation times for chunk window ")

    parser.add_argument('--partition_size', type=int, default=15,
                        help="EXPERIMENTAL: Maximum variants in per group size")

    parser.add_argument('--minimum_allele_gap', type=int, default=0.15,
                        help="EXPERIMENTAL: Minimum allele gap filtering candidate path generation")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, will choose candidate +/- 1 or +/- 2. Use together with gen4Training. default: %(default)s")

    # options for internal process control
    ## Subsample ratio tag for sub-sampled BAM file
    parser.add_argument('--subsample_ratio', type=int, default=1000,
                        help=SUPPRESS)

    ## Test in specific candidate site. Only use for analysis
    parser.add_argument('--test_pos', type=str2bool, default=True,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    args = parser.parse_args()
    UnifyRepresentation(args)


if __name__ == "__main__":
    main()

# /autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/pypy3 /autofs/bal33/zxzheng/ru/ru.py UnifyRepresentation