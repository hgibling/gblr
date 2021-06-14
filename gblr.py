#! /usr/bin/env python

### import libraries
import argparse
import edlib
import math
import pandas as pd
import pysam
import os
import sys

### define functions
def get_sequences_from_fasta(file_name):
    out = list()
    file_handle = pysam.FastxFile(file_name)
    for sequence in file_handle:
        out.append(sequence)
    return out

def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, 'N') for base in reversed(sequence))
    return reverse_complement

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alleles', type=str, required=True, help='fasta file of allele sequences with flanking sequences')
parser.add_argument('-r', '--reads', type=str, required=True, help='fasta/q file of sequencing reads')
parser.add_argument('-f', '--flank-length', type=int, default=10000, help='length of sequences flanking alleles')
parser.add_argument('-e', '--max-edit-distance', type=int, default=20, help='maximum edit distance allowed for read corrections')
parser.add_argument('-d', '--debug', type=int, default=0, help='toggle printing of content for debugging')
args = parser.parse_args()

### get allele sequences, reads, and flanking sequences length
alleles = get_sequences_from_fasta(args.alleles)
reads = pysam.FastxFile(args.reads)
flank_length = args.flank_length

### define variables
num_quality_reads = 0
num_flank_reads = 0
num_bad_reads = 0
allele_counts = {}
edit_distances = []

### align each read against all alleles
for read in reads:
    best_distance = args.max_edit_distance
    best_allele = []
    start = []
    end = []
    
    ### check alignment for forward and reverse reads
    for strand_idx, strand_sequence in enumerate([read.sequence, reverse_complement(read.sequence)]):
        ### check alignment to each allele
        for allele_idx, allele in enumerate(alleles):
            allele_length_no_end_flank = len(allele.sequence) - flank_length
            result = edlib.align(strand_sequence, allele.sequence, mode = "HW", task = "path")
            ### update variables if edit distance is acceptable or improved
            # ignore alignments that start after the last 50bp or end before the first 50bp of the variable region of interest
            if (result['locations'][0][0] > (allele_length_no_end_flank - 50)) | (result['locations'][0][1] < (flank_length + 50))
                num_flank_reads += 1
                continue
            if result['editDistance'] < best_distance:
                best_allele = [allele.name]
                best_distance = result['editDistance']
                start = [result['locations'][0][0]]
                end = [result['locations'][0][1]]
            elif result['editDistance'] == best_distance:
                best_allele.append(allele.name)
                start.append(result['locations'][0][0])
                end.append(result['locations'][0][1])
            else:
                num_bad_reads += 1

    ### consider only reads with acceptable edit distances to an allele
    ### case for single best allele
    if len(best_allele) ==  1:
        num_quality_reads += 1
        edit_distances.append(best_distance)
        allele_counts[best_allele] = allele_counts.get(best_allele, 0) + 1
        
    ### case where there are 2 or more equally best alleles
    elif len(best_allele) > 1:
        num_quality_reads += 1  
        edit_distances.append(best_distance)

        # iterate through all best alleles and count proportions
        for i, allele in enumerate(best_allele):
            allele_counts[best_allele[i]] = allele_counts.get(best_allele[i], 0) + 1/len(best_allele)

### print useful information
print("Edit distance used: %d" % args.max_edit_distance, file=sys.stderr)
print("Number of alleles tested: %d" % )
print("Number of quality reads: %d" % num_quality_reads, file=sys.stderr)
print("Number of low quality reads: %d" % num_bad_reads, file=sys.stderr)
print("Number of reads aligning predominantly to flank sequences: %d" % num_flank_reads, file=sys.stderr)

### print resulting counts for each allele
for allele, count in allele_counts.items():
    print(allele, '\t', count)
allele_count == 0:
