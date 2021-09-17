#! /usr/bin/env python

### import libraries
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline     # requires MAFFT to be installed and in PATH
from collections import defaultdict
from io import StringIO

import argparse
import edlib
import math
import numpy as np
import os
import pandas as pd
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

# get multiple sequence alignment
def get_MSA(allele, allele_reads_list, all_subset_reads):
    # get read sequences for each allele
    temp_genotype_reads_dict = dict.fromkeys(allele_reads_list, 0)
    for read in temp_genotype_reads_dict.keys():
        temp_genotype_reads_dict[read] = all_subset_reads[read]

    # write reads to temporary fasta
    fasta_name = "-".join([allele, "reads-temp.fa"])
    outfile = open(fasta_name, "w")
    for read, sequence in temp_genotype_reads_dict.items():
        outfile.write(">" + read + "\n")
        outfile.write(sequence + "\n")
    outfile.close

    # run mafft to get multiple sequence alignment
    mafft_command = MafftCommandline(input = fasta_name, localpair=True, maxiterate = 1000)
    stdout, stderr = mafft_command()

    # with open("aligned.fasta", "w") as handle:
    #     handle.write(stdout)
    # reads_MSA = AlignIO.read("aligned.fasta", "fasta")

    reads_MSA = AlignIO.read(StringIO(stdout), "fasta")
    return(reads_MSA)

# get consensus sequence (code modified from: https://stackoverflow.com/questions/38586800/python-multiple-consensus-sequences)
def get_consensus(reads_MSA):
    IUPAC_ambiguous = {'AG': 'R', 'CT': 'Y', 'CG': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M', 'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'ACGT': 'N'}

    alignment_length = reads_MSA.get_alignment_length()
    profile = {'A': [0]*alignment_length, 'C': [0]*alignment_length, 'G': [0]*alignment_length, 'T': [0]*alignment_length, '-': [0]*alignment_length}

    # count nucleotide occurances at each location
    for record in reads_MSA:
        for i, nuc in enumerate(record.seq.upper()):
            profile[nuc][i] += 1
    
    # determine consensus sequence
    consensus = ""
    for i in range(alignment_length):
        max_count = 0
        max_nuc = ""
        for nuc in "ACGT-":
            if profile[nuc][i] > max_count:
                max_count = profile[nuc][i]
                max_nuc = nuc
            elif profile[nuc][i] == max_count:
                max_nuc = "".join([max_nuc, nuc])
        if len(max_nuc) > 1:
            if '-' in max_nuc:
                max_nuc = max_nuc.replace("-", "")
                if len(max_nuc) == 1:
                    max_nuc = max_nuc.lower()
                else: 
                    max_nuc = IUPAC_ambiguous[max_nuc].lower()
            else:
                max_nuc = IUPAC_ambiguous[max_nuc]
        consensus += max_nuc
    return(consensus.replace("-", ""))

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--alleles', type=str, required=True, help='fasta file of sequences for alleles or region of interest, including flanking sequences')
parser.add_argument('-r', '--reads', type=str, required=True, help='bam file of aligned sequencing reads (or fastx if --quick-count is specified)')
parser.add_argument('-R', '--region', type=str, default="5:23526782,23527873", help='position of one region of interest with a chromosome name that matches the bam provided for --reads (ex: chr1:100000-200000)')
parser.add_argument('-l', '--flank-length', type=int, default=10000, help='length of sequences flanking alleles')
parser.add_argument('-t', '--flank-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the flanking regions')
parser.add_argument('-e', '--error-rate', type=float, default=0.001, help='estimate of the sequencing error rate')
parser.add_argument('-d', '--diploid', action='store_true', help='get diploid genotype scores instead of haploid (cannot be used with --quick-count)')
parser.add_argument('-q', '--quick-count', action='store_true', help='get counts of reads that align best to alleles instead of scores')
parser.add_argument('-m', '--max-mismatch', type=float, default=0.05, help='for quick count: maximum proportion of a read that can be mismatched/indels relative to an allele')
parser.add_argument('-T', '--alignment-tolerance', type=int, default=50, help='for quick count: minimum number of bases to which a read must align in the variable region of interest')
parser.add_argument('-D', '--delimiter', type=str, default='\t', help='delimiter to use for results output')
parser.add_argument('-v', '--verbose', action='store_true', help='print table of edit distances to stderr')
parser.add_argument('-A', '--alignments', type=str, help='print alignments of specified alleles (separated by commas; ex: C,D,L18) plus best allele to stderr')
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
quality_reads = set()
flank_reads = set()
region_of_interest_reads = set()
bad_reads = set()
allele_names = list(alleles.keys())
all_allele_lengths = dict.fromkeys(allele_names)
edit_distances = []

### make sure specified flank length not longer than any allele sequences
for name, sequence in alleles.items():
    if len(sequence) <= args.flank_length * 2:
        exit("ERROR: flank-length argument too long for allele sequence %s" % name)
    else:
        all_allele_lengths[name] = len(sequence)

### for generating quick counts
if args.quick_count:
    allele_counts = defaultdict(int)

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
    all_subset_reads = {}
    region = re.split("[:,-]", args.region)
    region = [region[0], int(region[1]), int(region[2])]
    all_alignments = {}
    if args.alignments != None:
        alignment_alleles = args.alignments.split(",")

    ### iterate over reads that overlap with region of interest
    for read in reads.fetch(region[0], region[1], region[2]):

        ### filter reads to consider only those that touch both flanks
        if read.reference_start < (region[1] - args.flank_tolerance) and read.reference_end > (region[2] + args.flank_tolerance):
            read_distance_dict = dict.fromkeys(allele_names)
            read_alignment_dict = dict.fromkeys(allele_names)
                
            ### subset read to just the region of interest
            chop_sites = subset_positions(read.cigarstring, read.reference_start, read.reference_end, region[1], region[2]) 
            read_subset = read.query_alignment_sequence[chop_sites[0] : -chop_sites[1]]
            all_subset_reads[read.query_name] = read_subset
            
            ### check alignment to each allele
            for allele_name, allele_sequence in alleles.items():
                best_distance = math.inf
            
                ### check alignment for forward and reverse reads
                for strand_idx, strand_sequence in enumerate([read_subset, reverse_complement(read_subset)]):
                    subset_alignment = edlib.align(strand_sequence, allele_sequence[args.flank_length : -args.flank_length], mode = "NW", task = "path")
                    region_of_interest_reads.add(read.query_name) 

                    ### check edit distance and store lowest edit distance between forward and reverse read sequences
                    if subset_alignment['editDistance'] <= best_distance:
                        read_distance_dict[allele_name] = subset_alignment['editDistance']
                        best_distance = subset_alignment['editDistance']
                        read_alignment_dict[allele_name] = edlib.getNiceAlignment(subset_alignment, strand_sequence, allele_sequence[args.flank_length : -args.flank_length])
            
            ### store read edit distances for each allele
            all_edit_distances[read.query_name] = read_distance_dict
            all_alignments[read.query_name] = read_alignment_dict

    ### get table of edit distances         # TODO: deal with null values
    allele_edit_distances = pd.DataFrame.from_dict(all_edit_distances, orient='index')

    ### for each read, get allele with best alignment
    if args.alignments != None:
        best_alleles = allele_edit_distances.idxmin(axis=1)

    if args.verbose:
        ### get sum of edit distances to print
        allele_edit_distances.loc['Sum'] = allele_edit_distances.sum()
        print(allele_edit_distances.to_string(), file=sys.stderr)
        ### remove sum row for further analysis
        allele_edit_distances = allele_edit_distances[:-1]

    ### get genotype edit distances if doing diploid calling
    if args.diploid:

        ### get all possible genotypes
        genotype_names = get_genotype_names(allele_names)
        genotype_edit_distances = pd.DataFrame(index=allele_edit_distances.index, columns=genotype_names)

        for g in genotype_names:
            split_alleles = g.split('/')
            genotype_edit_distances[g] = np.logaddexp((allele_edit_distances[split_alleles[0]] * np.log(args.error_rate) - np.log(2)), (allele_edit_distances[split_alleles[1]] * np.log(args.error_rate) - np.log(2)))

        all_scores = genotype_edit_distances.sum().sort_values(ascending=False)

        ### from top genotype, find out which read is most likely from which allele
        top_genotype_split = all_scores.index[0].split('/')
        reads_best_allele = allele_edit_distances[top_genotype_split].idxmin(axis=1)

        top_genotype_subset_reads = {}
        for a in top_genotype_split:
            top_genotype_subset_reads[a] = list(reads_best_allele[reads_best_allele==a].index)

        ### get consensus sequences of the reads for each allele in the top genotype
        novel_counter = 0
        for allele in top_genotype_subset_reads.keys():
            allele_MSA = get_MSA(allele, top_genotype_subset_reads[allele], all_subset_reads)
            allele_consensus = get_consensus(allele_MSA)
            allele_subsequence = alleles[allele][args.flank_length:-args.flank_length]
            if allele_consensus != allele_subsequence:
                novel_counter += 1
                print("Read consensus sequence for allele %s in top-scoring genotype does not match allele sequence: sample likely a has a novel haplotype" % allele, file=sys.stderr)
                #print("Read consensus:   %s" % allele_consensus, file=sys.stderr)
                #print("Allele consensus: %s" % allele_subsequence, file=sys.stderr)
        if novel_counter == 1:
            print("Heterozygous-Novel", "NA", sep=args.delimiter)
        elif novel_counter == 2:
            print("Homozygous-Novel", "NA", sep=args.delimiter)

    else:   # haploid calling
        all_scores = allele_edit_distances.sum().sort_values()
    
    ### print results (allele or genotype name, score)
    for name, score in all_scores.items():
        print(name, score, sep=args.delimiter)

    ### for each read, print alignments for best allele and specified alleles
    if args.alignments != None:
        for read, dictionary in all_edit_distances.items():
            print("Read %s: best alelle was %s (ED=%d)" % (read, best_alleles[read], dictionary.get(best_alleles[read])), file=sys.stderr)
            print("\n".join(all_alignments[read][best_alleles[read]].values()), file=sys.stderr)
            print("\n", file=sys.stderr)
            for a in alignment_alleles:
                print("Compare to alginment to allele %s (ED=%d):" % (a, dictionary[a]), file=sys.stderr)
                print("\n".join(all_alignments[read][a].values()), file=sys.stderr)
                print("\n", file=sys.stderr)
            print("---\n", file=sys.stderr)

### print useful information to stderr
# print("Max proportion of read that is mismatches/indels: %d" % args.max_mismatch, file=sys.stderr)
# print("Edit distance stats: Min: %d, Median: %d, Max: %d" % (ED_min, ED_mean, ED_max), file=sys.stderr)
# print("Number of alleles tested: %d" % len(alleles), file=sys.stderr)
# print("Number of quality reads: %d" % len(quality_reads), file=sys.stderr)
# print("Number of low quality reads: %d" % len(bad_reads), file=sys.stderr)
# print("Number of reads aligning predominantly to flank sequences: %d" % len(flank_reads), file=sys.stderr)