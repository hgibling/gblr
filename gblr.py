#! /usr/bin/env python

### import libraries
from Bio import AlignIO
from collections import defaultdict

import argparse
import edlib
import math
import numpy as np
import os
import pandas as pd
import pyabpoa as pa
import pysam
import re
import sys

### define functions
# read in fastx file and save name and sequence from each entry
def get_sequences_from_fasta(file_name):
    out = {}
    file_handle = pysam.FastxFile(file_name)
    for fasta in file_handle:
        out[fasta.name] = fasta.sequence
    return out

# get reverse complement of a sequence
def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, 'N') for base in reversed(sequence))
    return reverse_complement

# get all possible genotypes
def get_genotype_names(allele_names):
    genotype_names = []
    for i in allele_names:
        for j in allele_names:
            if i <= j:
                genotype_names.append(i + '/' + j)
    return genotype_names

# figure how much to adjust left and right chops for subsetting reads when they have indels before/after the region of interest
def subset_positions(read_cigar, ref_start, ref_end, region_start, region_end):
    left_chop = region_start - ref_start - 1
    right_chop = ref_end - region_end
    split_cigar = re.findall('[0-9]*[A-Z=]', read_cigar)
    left_count, left_indel, right_count, right_indel = 0, 0, 0, 0

    # adjust left_chop to account for indels in read before the region of interest
    for chunk in split_cigar: 
        if left_count >= left_chop:
            break   # start of the region of interest has been reached
        if chunk.endswith(('M', 'X', '=')):
            left_count += int(chunk[:chunk.find('MX=')])
        elif chunk.endswith('D'):
            left_indel -= int(chunk[:chunk.find('D')])
            left_count -= int(chunk[:chunk.find('D')])
        elif chunk.endswith('I'):
            left_indel += int(chunk[:chunk.find('I')])
            left_count += int(chunk[:chunk.find('I')])
    left_chop += left_indel
    
    # adjust right_chop to account for indels in read after the region of interest
    for chunk in split_cigar[::-1]: 
        if right_count >= right_chop:
            break   # start of the region of interest has been reached
        if chunk.endswith(('M', 'X', '=')):
            right_count += int(chunk[:chunk.find('MX=')])
        elif chunk.endswith('D'):
            right_indel -= int(chunk[:chunk.find('D')])
            right_count -= int(chunk[:chunk.find('D')])
        elif chunk.endswith('I'):
            right_indel += int(chunk[:chunk.find('I')])
            right_count += int(chunk[:chunk.find('I')])
    right_chop += right_indel

    # return values for subsetting
    return [left_chop, right_chop]

# modify cigar sequence so indels only count as an edit distance of 1
def modify_edit_distance(cigar):
    split_cigar = re.findall('[0-9]*[A-Z=]', cigar)
    edit_distance = 0
    for chunk in split_cigar:
        if chunk.endswith('X'):
            edit_distance += int(chunk[:chunk.find('X')])
        elif chunk.endswith(('D', 'I')):
            edit_distance += 1
    return(edit_distance)

# get multiple sequence alignment with ABPOA
def get_MSA(allele_reads_list, all_subset_reads):
    # define aligner parameters
    aligner = pa.msa_aligner()
    # get reads of interest
    reads = [all_subset_reads[key] for key in allele_reads_list]
    # align all reads
    reads_MSA = aligner.msa(reads, out_cons=True, out_msa=True)
    # return dict of list of seqs in MSA format, as well consensus sequence
    # dict{msa: [reads], consensus: [consensusseq]}
    return({'MSA':reads_MSA.msa_seq, 'consensus':reads_MSA.cons_seq})

# get consensus alignment with ABPOA
def get_abpoa_consensus(all_reads):
    # define aligner parameters
    aligner = pa.msa_aligner()
    # align all reads
    reads_consensus = aligner.msa(list(all_subset_reads.values()), out_cons=True, out_msa=False, max_n_cons=2, min_freq=0.35)
    # return list of seqs in MSA format, as well as consensus sequence
    # list[consensusseq]
    return(reads_consensus.cons_seq)
    
# all IUPAC ambiguous nucleotides
IUPAC_ambiguous_to_nucleotides = {'R':'AG', 'Y':'CT', 'S':'CG', 'W':'AT', 'K':'GT', 'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG', 'N':'ACGT'}

# get consensus sequence (code modified from: https://stackoverflow.com/questions/38586800/python-multiple-consensus-sequences)
def get_consensus(reads_MSA, threshold = 0.35, ambiguous = IUPAC_ambiguous_to_nucleotides):
    IUPAC_nucs_to_ambiguous = dict((nuc, amb) for amb, nuc in ambiguous.items())

    # alignment_length = reads_MSA.get_alignment_length()
    alignment_length = len(reads_MSA[1])
    profile = pd.DataFrame({'A': [0]*alignment_length, 'C': [0]*alignment_length, 'G': [0]*alignment_length, 'T': [0]*alignment_length, '-': [0]*alignment_length})

    # count nucleotide occurances at each location
    # for record in reads_MSA:
    #     for i, nuc in enumerate(record.seq.upper()):
    #         profile[nuc][i] += 1
    for seq in reads_MSA:
        for i, nuc in enumerate(seq):
            profile[nuc][i] += 1
    profile = profile.transpose()
    
    # determine consensus sequence
    consensus = ""
    for i in range(alignment_length):
        # get all nucleotides that have a frequency of at least defined threshold
        max_nuc = "".join(list(profile[profile[i] > (profile[i].sum()*threshold)].index))
        if len(max_nuc) > 1:
            if '-' in max_nuc:
                max_nuc = max_nuc.replace("-", "")
                if len(max_nuc) == 1:
                    # use lowercase letters to notate positions where a gap occurs ~as frequently as the letter
                    max_nuc = max_nuc.lower()
                else: 
                    max_nuc = IUPAC_nucs_to_ambiguous[max_nuc].lower()
            else:
                max_nuc = IUPAC_nucs_to_ambiguous[max_nuc]
        consensus += max_nuc
    return(consensus.replace("-", ""))

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alleles', type=str, required=True, help='fasta file of sequences for alleles or region of interest, including flanking sequences')
parser.add_argument('-r', '--reads', type=str, required=True, help='bam file of aligned sequencing reads (or fastx if --quick-count is specified)')
parser.add_argument('-R', '--region', type=str, default="5:23526673,23527764", help='position of the region of interest (ex: chr1:100000-200000)')
parser.add_argument('-l', '--flank-length', type=int, default=10000, help='length of sequences flanking alleles')
parser.add_argument('-t', '--flank-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the flanking regions')
parser.add_argument('-e', '--error-rate', type=float, default=0.01, help='estimate of the sequencing error rate')
parser.add_argument('-d', '--diploid', action='store_true', help='get diploid genotype scores instead of haploid (cannot be used with --quick-count)')
parser.add_argument('-s', '--scoring-model', type=str, help='scoring model to use ("e" or "1ee")')
parser.add_argument('-E', '--ED-model', type=str, help='edit distance model to use ("allIndel" or "1Indel")')
parser.add_argument('-N', '--print-top-N-genos', type=int, default=0, help='print likelihoods of only the top N genotypes (default: print all')
parser.add_argument('-v', '--verbose', action='store_true', help='print table of edit distances to stderr')
parser.add_argument('-V', '--verbose-reads', action='store_true', help='print list of reads that best align to each allele in top genotype')
parser.add_argument('-c', '--consensus_sequence', action='store_true', help='print consensus sequences to output-name.consensus.fa')
# parser.add_argument('-C', '--consensus_alignment', action='store_true', help='print alignment of read consensus sequence and alleles from top genotype to stderr') # TODO
parser.add_argument('-q', '--quick-count', action='store_true', help='get counts of reads that align best to alleles instead of scores')
parser.add_argument('-m', '--max-mismatch', type=float, default=0.05, help='for quick count: maximum proportion of a read that can be mismatched/indels relative to an allele')
parser.add_argument('-T', '--alignment-tolerance', type=int, default=50, help='for quick count: minimum number of bases to which a read must align in the variable region of interest')
parser.add_argument('-D', '--delimiter', type=str, default='\t', help='delimiter to use for results output')
parser.add_argument('-o', '--output-name', type=str, required=True, help='name of file to save scores/calls')
args = parser.parse_args()

### check arguments
if args.flank_length < args.flank_tolerance:
    exit("ERROR: flank length must be longer than flank tolerance. Check parameters.")

if args.max_mismatch < 0 or args.max_mismatch > 1:
    exit("ERROR: max-mistmatch argument must be a proportion (between 0 and 1). Check parameters.")

if args.diploid and args.quick_count:
    exit("ERROR: diploid calling cannot be done with the quick-count methid (haploid only). Check parameters.")

### get allele and read sequences
alleles = get_sequences_from_fasta(args.alleles)
if args.quick_count:
    with open(args.reads, 'r') as temp:
        try:
            temp.readline()   
            reads = pysam.FastxFile(args.reads)
        except:
            exit("ERROR: issue with the reads file. Is it a proper fasta or fastq file?")
else:
    try:
        reads = pysam.AlignmentFile(args.reads, 'rb')
    except:
        exit("ERROR: issue with the reads file. Is it a proper bam file?")

### define variables
region_of_interest_reads = set()
allele_names = list(alleles.keys())
all_allele_lengths = dict.fromkeys(allele_names)
edit_distances = []
N_geno = 0

### make sure specified flank length not longer than any allele sequences
for name, sequence in alleles.items():
    if len(sequence) <= args.flank_length * 2:
        exit("ERROR: flank-length argument too long for allele sequence %s" % name)
    else:
        all_allele_lengths[name] = len(sequence)

### set up output file
results_file = open(args.output_name, "w")

### for generating quick counts
if args.quick_count:
    quality_reads = set()
    flank_reads = set()
    bad_reads = set()
    allele_counts = defaultdict(int)

    ### align each read against all alleles
    for read in reads:
        best_distance = args.max_mismatch * len(read.sequence)
        best_allele = []
        
        ### check alignment for forward and reverse reads
        for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
            ### check alignment to each allele
            for allele_name, allele_sequence in alleles.items():
                result = edlib.align(strand_sequence, allele_sequence, mode = "NW", task = "path")

                ### ignore alignments that do not meet the minimum number of bases a read must align to in the variable region of interest
                # check if alignment starts after the 3' minimum alignment threshold, or ends before the 5' minimum threshold
                if ((result['locations'][0][0] > (all_allele_lengths[allele_name] - args.flank_length - args.alignment_tolerance)) or (result['locations'][0][1] < (args.flank_length + args.alignment_tolerance))):
                    flank_reads.add(read.name)
                    continue

                ### consider alignments that meet the maximum edit distance threshold
                # update variables if edit distance is under acceptable threshold or is improved from a previous acceptable alignment
                if result['editDistance'] < best_distance:
                    best_allele = [allele_name]
                    best_distance = result['editDistance']
                    quality_reads.add(read.name)
                # add to list of best alleles if edit distance is the same as a previously determined acceptable alignment
                elif result['editDistance'] == best_distance:
                    best_allele.append(allele_name)
                    quality_reads.add(read.name) 

        ### get metrics for read classes
        # if alignment to one allele was quality but to the flank sequences in another allele, only consider quality alignment 
        if read.name in flank_reads and quality_reads:
            flank_reads.remove(read.name)
        # if read didn't align to any part of the allele sequence under the edit distance threshold, consider it low quality
        if (read.name not in flank_reads) and (read.name not in quality_reads):
            bad_reads.add(read.name)

        ### consider only reads with acceptable edit distances to an allele and add count to best allele(s)
        # single best allele
        if len(best_allele) ==  1:
            edit_distances.append(best_distance)
            allele_counts[best_allele[0]] += 1    
                
        # two or more equally best alleles
        elif len(best_allele) > 1:
            edit_distances.append(best_distance)
            # iterate through all best alleles and count proportions
            for i, allele in enumerate(best_allele):
                allele_counts[best_allele[i]] += 1/len(best_allele)
        
        else:
            exit("ERROR: no reads aligned to region of interest")

    ### print results (allele, count, proportion)
    for allele, count in sorted(allele_counts.items(), key=lambda x: x[1], reverse=True):
        print(allele, count, count/len(quality_reads), sep=args.delimiter, file=results_file)
        N_geno += 1
        if args.print_top_N_genos > 0 and N_geno == args.print_top_N_genos:
            break

### for generating scores:
else:
    all_edit_distances = {}
    all_subset_reads = {}
    region = re.split("[:,-]", args.region)
    region = [region[0], int(region[1]), int(region[2])]
    all_alignments = {}

    ### check chromosome naming convention and update region if needed
    if ("chr" in region[0]) != ("chr" in reads.references[0]):
        if "chr" in reads.references[0]:
            region[0] = "chr" + region[0]
        else:
            region[0] = region[0].replace("chr", "")
    
    ### iterate over reads that overlap with region of interest
    for read in reads.fetch(region[0], region[1], region[2]):

        ### filter reads to consider only those that touch both flanks
        if read.reference_start < (region[1] - args.flank_tolerance) and read.reference_end > (region[2] + args.flank_tolerance):
            read_distance_dict = dict.fromkeys(allele_names)
                
            ### subset read to just the region of interest
            chop_sites = subset_positions(read.cigarstring, read.reference_start, read.reference_end, region[1], region[2])
            try:
                read_subset = read.query_alignment_sequence[chop_sites[0] : -chop_sites[1]]
            except TypeError:   # to ignore reads without sequences
                continue
            all_subset_reads[read.query_name] = read_subset
            
            ### check alignment to each allele
            for allele_name, allele_sequence in alleles.items():
                subset_alignment = edlib.align(read_subset, allele_sequence[args.flank_length : -args.flank_length], mode = "NW", task = "path")
                region_of_interest_reads.add(read.query_name)
                if (args.ED_model == "1Indel"):
                    read_distance_dict[allele_name] = modify_edit_distance(subset_alignment['cigar'])
                elif (args.ED_model == "allIndel"):
                    read_distance_dict[allele_name] = subset_alignment['editDistance']
            
            ### store read edit distances for each allele
            all_edit_distances[read.query_name] = read_distance_dict
    
    ### print number of reads used for analysis
    print("Number of reads that fully span region of interest: %d" % (len(region_of_interest_reads)), file=sys.stderr)

    ### get table of edit distances
    # dataframe[reads,alleles: edit distance]         # TODO: deal with null values
    allele_edit_distances = pd.DataFrame.from_dict(all_edit_distances, orient='index')

    if args.verbose:
        ### get sum of edit distances to print
        allele_edit_distances.loc['Sum'] = allele_edit_distances.sum()
        print(allele_edit_distances.to_string(), file=sys.stderr)
        ### remove sum row for further analysis
        allele_edit_distances = allele_edit_distances[:-1]

    ### get genotype edit distances if doing diploid calling
    if args.diploid:

        ### get read subset lengths and add to edit distance dict
        all_subset_reads_dict = pd.DataFrame.from_dict(all_subset_reads, orient='index', columns=['ReadLength'])
        all_subset_reads_dict.ReadLength = all_subset_reads_dict.ReadLength.str.len()
        allele_edit_distances['ReadLength'] = all_subset_reads_dict

        ### get all possible genotypes
        genotype_names = get_genotype_names(allele_names)
        genotype_edit_distances = pd.DataFrame(index=allele_edit_distances.index, columns=genotype_names)

        ### get genotype likelihoods
        # dataframe[reads,genos: likelihoods]
        if args.scoring_model == "e":
            for g in genotype_names:
                split_alleles = g.split('/')
                genotype_edit_distances[g] = np.logaddexp(((allele_edit_distances[split_alleles[0]] * np.log(args.error_rate)) - np.log(2)), ((allele_edit_distances[split_alleles[1]] * np.log(args.error_rate)) - np.log(2)))
        elif args.scoring_model == "1ee":
            for g in genotype_names:
                split_alleles = g.split('/')
                genotype_edit_distances[g] = np.logaddexp((((allele_edit_distances['ReadLength'] - allele_edit_distances[split_alleles[0]]) * np.log(1 - args.error_rate)) + (allele_edit_distances[split_alleles[0]] + np.log(args.error_rate)) - np.log(2)), (((allele_edit_distances['ReadLength'] - allele_edit_distances[split_alleles[1]]) * np.log(1 - args.error_rate)) + (allele_edit_distances[split_alleles[1]] + np.log(args.error_rate)) - np.log(2)))

        ### get overall likelihood for each genotype
        # series[genos: likelihood]
        all_scores = genotype_edit_distances.sum().sort_values(ascending=False)
        top_genotype_split = list(set(all_scores.index[0].split('/')))

        ### get abpoa consensus
        consensus_seqs = get_abpoa_consensus(all_subset_reads)

        ### print stats
        print("Top genotype is %s" % (all_scores.index[0]), file=sys.stderr)
        print("Number of consensus sequences is %d" % (len(consensus_seqs)), file=sys.stderr)
        if len(consensus_seqs) != len(top_genotype_split):
            print("Top genotype zygosity does not match consensus sequence(s)", file=sys.stderr)
        if args.consensus_sequence:
            # print fasta for consensus sequence(s)
            consensus_file = open(args.output_name + ".consensus.fa", "a")
            print(">consensus sequence", file=consensus_file)
            print(consensus_seqs[0], file=consensus_file)
            if len(consensus_seqs) > 1:
                print(">second consensus sequence", file=consensus_file)
                print(consensus_seqs[1], file=consensus_file)
            consensus_file.close

        # 9 possible outcomes:
        # Num cons seqs | Num alleles | Num matches |
        # --------------|-------------|-------------|
        #             1 |           1 |          1  |   ex. A/A
        #             1 |           1 |          0  |   ex. NovelA/NovelA
        #             2 |           1 |          1  |   ex. NovelA/A
        #             2 |           1 |          0  |   ex. NovelA.a/NovelA.b
        #             1 |           2 |          1  |   ex. A/PotentialFalseB
        #             1 |           2 |          0  |   ex. PotentialFalseOrNovelA/PotentialFalseOrNovelB
        #             2 |           2 |          2  |   ex. A/B
        #             2 |           2 |          1  |   ex. NovelA/B
        #             2 |           2 |          0  |   ex. NovelA/NovelB

        novel_alleles = []
        known_alleles = []

        for allele in top_genotype_split:
            allele_subsequence = alleles[allele][args.flank_length:-args.flank_length]
            # check for novel haplotypes
            if consensus_seqs[0] == allele_subsequence:
                known_alleles.append(allele)
            elif len(consensus_seqs) == 2:
                if consensus_seqs[1] == allele_subsequence:
                    known_alleles.append(allele)
                else:
                    novel_alleles.append(allele)
            else: 
                novel_alleles.append(allele)

        # if there are any novel alleles detected, add them to the top of the likelihood results
        # NOTE: value of 1 is to ensure novel genotype stays at the top of the list--it is not a likelihood score            
        if len(novel_alleles) == 1:
            novel_name1 = "_".join(["Novel_Similar", novel_alleles[0]])
            if len(known_alleles) == 1:
                # only possible if top geno is het
                if len(consensus_seqs) == 1:
                    # unexpected result: give arbitrary value of 99
                    extra_name1 = "_".join(["Potential_False", novel_alleles[0]])
                    print("/".join([known_alleles[0], extra_name1]), "99", sep=args.delimiter, file=results_file)
                elif len(consensus_seqs) == 2:
                    # two base alleles, one novel and one known (ex. NovelA/B)
                    print("/".join([novel_name1, known_alleles[0]]), "1", sep=args.delimiter, file=results_file)
            elif len(known_alleles) == 0:
                # only possible if top geno is hom
                if len(consensus_seqs) == 1:
                    # one base allele, both same novel (ex. NovelA/NovelA)
                    print("/".join([novel_name1, novel_name1]), "1", sep=args.delimiter, file=results_file)
                elif len(consensus_seqs) == 2:
                    # one base allele, two different novel alleles (ex. NovelA.a/NovelA.b)
                    novel_name2 = "_".join(["Different_Novel_Similar", novel_alleles[0]])
                    print("/".join([novel_name1, novel_name2]), "1", sep=args.delimiter, file=results_file)
        elif len(novel_alleles) == 2:
            # only possible if top geno is het
            if len(consensus_seqs) == 1:
                # unexpected result: give arbitrary value of 99
                extra_name1 = "_".join(["Potential_False_or_Novel", novel_alleles[0]])
                extra_name2 = "_".join(["Potential_False_or_Novel", novel_alleles[1]])
                print("/".join([extra_name1, extra_name2]), "99", sep=args.delimiter, file=results_file)
            elif len(consensus_seqs) == 2:
                # two base alleles, both novel (ex. NovelA/NovelB)
                novel_name1 = "_".join(["Novel_Similar", novel_alleles[0]])
                novel_name2 = "_".join(["Novel_Similar", novel_alleles[1]])
                print("/".join([novel_name1, novel_name2]), "1", sep=args.delimiter, file=results_file)  
        elif len(novel_alleles) == 0 and len(known_alleles) == 1 and len(consensus_seqs) == 2: 
            # only possible if top geno is hom
            # one base allele, one novel and one known (ex. NovelA/A)
            novel_name1 = "_".join(["Novel_Similar", known_alleles[0]])
            print("/".join([novel_name1, known_alleles[0]]), "1", sep=args.delimiter, file=results_file)
        # elif len(known_alleles) == 1 and len(consensus_seqs) == 1: A/A printed below
        # elif len(known_alleles) == 2: A/B printed below

    else:   # haploid calling
        all_scores = allele_edit_distances.sum().sort_values()
        # TODO: check for novel alleles
    
    ### print results (allele or genotype name, score)
    for name, score in all_scores.items():
        print(name, score, sep=args.delimiter, file=results_file)
        N_geno += 1
        if (args.print_top_N_genos > 0) and (N_geno == args.print_top_N_genos):
            break

results_file.close
