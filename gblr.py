#! /usr/bin/env python

### import libraries
from collections import defaultdict
import argparse
import edlib
import math
import numpy as np
import os
import pandas as pd
import pysam
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

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alleles', type=str, required=True, help='fasta file of sequences for alleles or region of interest, including flanking sequences')
parser.add_argument('-r', '--reads', type=str, required=True, help='fasta/q file of sequencing reads')
parser.add_argument('-f', '--flank-length', type=int, default=10000, help='length of sequences flanking alleles')
parser.add_argument('-t', '--flank-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the flanking regions')
parser.add_argument('-d', '--diploid', action='store_true', help='get diploid genotype scores instead of haploid (cannot be used with --quick-count)')
parser.add_argument('-q', '--quick-count', action='store_true', help='get counts of reads that align best to alleles instead of scores')
parser.add_argument('-m', '--max-mismatch', type=float, default=0.05, help='for quick count: maximum proportion of a read that can be mismatched/indels relative to an allele')
parser.add_argument('-T', '--alignment-tolerance', type=int, default=50, help='for quick count: minimum number of bases to which a read must align in the variable region of interest')
parser.add_argument('-D', '--delimiter', type=str, default='\t', help='delimiter to use for results output')
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
reads = pysam.FastxFile(args.reads)

### define variables
quality_reads = set()
flank_reads = set()
bad_reads = set()
allele_names = list(alleles.keys())
all_allele_lengths = dict.fromkeys(allele_names)
allele_counts = defaultdict(int)
edit_distances = []

### make sure specified flank length not longer than any allele sequences
for name, sequence in alleles.items():
    if len(sequence) <= args.flank_length * 2:
        exit("ERROR: flank-length argument too long for allele sequence %s" % name)
    else:
        all_allele_lengths[name] = len(sequence)

### for generating quick counts
if args.quick_count:

    ### align each read against all alleles
    for read in reads:
        best_distance = args.max_mismatch * len(read.sequence)
        best_allele = []
        
        ### check alignment for forward and reverse reads
        for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
            ### check alignment to each allele
            for allele_name, allele_sequence in alleles.items():
                result = edlib.align(strand_sequence, allele_sequence, mode = "HW", task = "path")

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

    ### calculate stats for edit distances between the quality reads and the best alleles
    ED_max = max(edit_distances)
    ED_min = min(edit_distances)
    ED_mean = sum(edit_distances) / len(edit_distances)

    ### print results (allele, count, proportion)
    for allele, count in sorted(allele_counts.items(), key=lambda x: x[1], reverse=True):
        print(allele, count, count/len(quality_reads), sep=args.delimiter)

### for generating scores:
else:
    all_edit_distances = {}
    all_cigars = {}

    ### align each read against all alleles
    for read in reads:
        read_distance_dict = dict.fromkeys(allele_names)
        read_cigar_dict = dict.fromkeys(allele_names)

        ### check alignment to each allele
        for allele_name, allele_sequence in alleles.items():
            best_distance = math.inf
        
            ### check alignment for forward and reverse reads
            for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
                result = edlib.align(strand_sequence, allele_sequence, mode = "HW", task = "path")


                ### check that read spans full region of interest and at least args.flank_tolerance into both flank sequences 
                ### if not, ignore
                if result['locations'][0][0] > (args.flank_length - args.flank_tolerance) or result['locations'][0][1] < (all_allele_lengths[allele_name] - args.flank_length + args.flank_tolerance):
                    flank_reads.add(read.name)
                    continue

                ##### subset read to the portion that aligns to the region of interest
                subset_start_position = args.flank_length - result['locations'][0][0]
                subset_end_position = args.flank_length - (all_allele_lengths[allele_name] - result['locations'][0][1]) + 1
                strand_subset = strand_sequence[subset_start_position : -end]

                ### realign read subset to allele
                subset_result = edlib.align(strand_subset, allele_sequence[args.flank_length : -args.flank_length], mode = "HW", task = "path")

                ### check edit distance and store lowest edit distance between forward and reverse read sequences
                if result['editDistance'] <= best_distance:
                    read_distance_dict[allele_name] = result['editDistance']
                    read_cigar_dict[allele_name] = result['cigar']
                    best_distance = result['editDistance']
            
        ### if acceptable alignment, store read edit proportions for each allele
        all_edit_distances[read.name] = read_distance_dict

    ### get table of edit distances         # TODO: deal with null values
    allele_edit_distances = pd.DataFrame.from_dict(all_edit_distances, orient='index')
    print(allele_edit_distances, file=sys.stderr)

    ### get genotype edit distances if doing diploid calling
    if args.diploid:

        ### get all possible genotypes
        genotype_names = get_genotype_names(allele_names)
        genotype_edit_distances = pd.DataFrame(index=allele_edit_distances.index, columns=genotype_names)

        for g in genotype_names:
            split_alleles = g.split('/')
            genotype_edit_distances[g] = ( allele_edit_distances[split_alleles[0]] / 2 ) + ( allele_edit_distances[split_alleles[1]] / 2 )

        print(genotype_edit_distances, file=sys.stderr)

        all_scores = genotype_edit_distances.sum().sort_values()

    else:

        all_scores = allele_edit_distances.sum().sort_values()
    
    ### print results (allele or genotype name, score)
    for name, likelihood in all_scores.items():
        print(name, likelihood, sep=args.delimiter)

### print useful information to stderr
# print("Max proportion of read that is mismatches/indels: %d" % args.max_mismatch, file=sys.stderr)
# print("Edit distance stats: Min: %d, Median: %d, Max: %d" % (ED_min, ED_mean, ED_max), file=sys.stderr)
# print("Number of alleles tested: %d" % len(alleles), file=sys.stderr)
# print("Number of quality reads: %d" % len(quality_reads), file=sys.stderr)
# print("Number of low quality reads: %d" % len(bad_reads), file=sys.stderr)
# print("Number of reads aligning predominantly to flank sequences: %d" % len(flank_reads), file=sys.stderr)