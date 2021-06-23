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

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alleles', type=str, required=True, help='fasta file of allele sequences with flanking sequences')
parser.add_argument('-r', '--reads', type=str, required=True, help='fasta/q file of sequencing reads')
parser.add_argument('-f', '--flank-length', type=int, default=10000, help='length of sequences flanking alleles')
parser.add_argument('-t', '--alignment-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the variable region of interest')
parser.add_argument('-m', '--max-mismatch', type=float, default=0.05, help='maximum proportion of a read that can be mismatched/indels relative to an allele')
parser.add_argument('-d', '--diploid', action='store_true', help='call diploid genotypes instead of haploid alleles')
parser.add_argument('-q', '--quick-count', action='store_true', help='get counts of reads that align best to alleles instead of full likelihood scoring')
parser.add_argument('-D', '--delimiter', type=str, default='\t', help='delimiter to use for output results')
args = parser.parse_args()

### check arguments
if args.max_mismatch < 0 or args.max_mismatch > 1:
    exit("ERROR: max-mistmatch argument must be a proportion (between 0 and 1). Check parameters.")

### get allele and read sequences
alleles = get_sequences_from_fasta(args.alleles)
reads = pysam.FastxFile(args.reads)

### define variables
quality_reads = set()
flank_reads = set()
bad_reads = set()
allele_names = list(alleles.keys())
allele_lengths = dict.fromkeys(allele_names)
allele_counts = defaultdict(int)
edit_distances = []

### make sure specified flank length not longer than any allele sequences
for name, sequence in alleles.items():
    allele_length_no_end_flank = len(sequence) - args.flank_length
    if allele_length_no_end_flank <= 0:
        exit("ERROR: flank-length argument too long for allele sequences provided (failed at allele %s)." % allele_name)
    else:
        allele_lengths[name] = allele_length_no_end_flank

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
                if ((result['locations'][0][0] > (allele_lengths[allele_name] - args.alignment_tolerance)) or (result['locations'][0][1] < (args.flank_length + args.alignment_tolerance))):
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

### for generating full likelihood scores:
else:
    mismatched_proportions = {}

    ### align each read against all alleles
    for read in reads:
        read_distance_dict = dict.fromkeys(allele_names)

        ### check alignment to each allele
        for allele_name, allele_sequence in alleles.items():
            best_distance = args.max_mismatch * len(read.sequence)
        
            ### check alignment for forward and reverse reads
            for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
                result = edlib.align(strand_sequence, allele_sequence, mode = "HW", task = "path")

                ### ignore alignments that do not meet the minimum number of bases a read must align to in the variable region of interest
                # check if alignment starts after the 3' minimum alignment threshold, or ends before the 5' minimum threshold
                if ((result['locations'][0][0] > (allele_lengths[allele_name] - args.alignment_tolerance)) or (result['locations'][0][1] < (args.flank_length + args.alignment_tolerance))):
                    flank_reads.add(read.name)
                    continue

                ### if acceptable edit distance, store lowest proportion of mismatched read bases between forward and reverse read sequences
                if result['editDistance'] <= best_distance:
                    read_distance_dict[allele_name] = result['editDistance'] / len(read.sequence)
                    best_distance = result['editDistance']
            
        ### if acceptable alignment, store read edit proportions for each allele
        mismatched_proportions[read.name] = read_distance_dict

    ### get log sums of values          # TODO: deal with null values
    mismatched_likelihoods = pd.DataFrame.from_dict(mismatched_proportions, orient='index')
    mismatched_likelihoods = np.log10(mismatched_likelihoods).sum().sort_values(ascending=False)
    
    ### print results (read, allele, edit distance proportion)
    for allele, likelihood in mismatched_likelihoods.items():
        print(allele, likelihood, sep=args.delimiter)



### print useful information to stderr
print("Max proportion of read that is mismatches/indels: %d" % args.max_mismatch, file=sys.stderr)
#print("Edit distance stats: Min: %d, Median: %d, Max: %d" % (ED_min, ED_mean, ED_max), file=sys.stderr)
print("Number of alleles tested: %d" % len(alleles), file=sys.stderr)
print("Number of quality reads: %d" % len(quality_reads), file=sys.stderr)
print("Number of low quality reads: %d" % len(bad_reads), file=sys.stderr)
print("Number of reads aligning predominantly to flank sequences: %d" % len(flank_reads), file=sys.stderr)