#! /usr/bin/env python

import argparse
import edlib
import math
import pandas as pd
import pysam
import os
import sys

def get_sequences_from_fasta(fn):
    out = list()
    fh = pysam.FastxFile(fn)
    for s in fh:
        out.append(s)
    return out

def reverse_complement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(s))
    return reverse_complement

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
read_flank_counts = {}
edit_distances = {}

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
                read_flank_counts += 1
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

    ### keep only reads with acceptable edit distances to an allele
    # case for single best allele
    if len(best_allele) ==  1:
        num_quality_reads += 1
        
        # reads that only align to right flank region
        if start[0] > max(allele_znfs.end):
            read_flanks.append("right")

        # # reads that fully or partially align to left flank region
        # elif start[0] < flank_length:
        #     read_flanks.append("left")
        #     # reads that extend alignment in znf region
        #     if end[0] >= min(allele_znfs.start):
        #         for z in allele_znfs.index:
        #             if end[0] >= allele_znfs.iloc[z].end:
        #                 read_znfs.append(allele_znfs.iloc[z].znfs)
        #             # if read ends alignment here, skip remaining znf iterations
        #             elif (end[0] > allele_znfs.iloc[z].start) & (end[0] < allele_znfs.iloc[z].end):
        #                 read_znfs.append(allele_znfs.iloc[z].znfs)
        #                 break

        # # reads that start alignment in znf region
        # elif start[0] >= flank_length:
        #     for z in allele_znfs.index:
        #         # skip over znfs that occur before the read starts aligning
        #         if start[0] > allele_znfs.iloc[z].end:
        #             continue
        #         if end[0] > allele_znfs.iloc[z].end:
        #             read_znfs.append(allele_znfs.iloc[z].znfs)
        #         # if read ends alignment here, skip remaining znf iterations
        #         elif (end[0] > allele_znfs.iloc[z].start) & (end[0] <= allele_znfs.iloc[z].end):
        #             read_znfs.append(allele_znfs.iloc[z].znfs)
        #             break
        #     # reads that end alignment in right flank region
        #     if end[0] > max(allele_znfs.end):
        #         read_flanks.append("right")

        # get znf and flank counts
        # for z in read_znfs:
        #     read_znf_counts[z] = read_znf_counts.get(z, 0) + 1
        # for f in read_flanks:
        #     read_flank_counts[f] = read_flank_counts.get(f, 0) + 1

        # if len(read_znfs) > 0:
        #     num_znfs_in_reads.append(len(read_znfs))

        # if args.debug:
        #print("%s, %s, %d, %d, %s" % (read.name, best_allele, start, end, read_znfs))

    # case where there are 2 or more equally best alleles
    elif len(best_allele) > 1:
        num_quality_reads += 1  
        all_best_read_znfs = []        

        # iterate through all best alleles
        for i in range(len(best_allele)):
            allele_znfs = znfs[znfs.allele == best_allele[i]].reset_index(drop=True)
            
            # reads that only align to right flank region
            if start[i] > max(allele_znfs.end):
                read_flanks.append("right")
            
            # # reads that fully or partially align to left flank region
            # elif start[i] < flank_length:
            #     read_flanks.append("left")
            #     # reads that extend alignment in znf region
            #     if end[i] >= min(allele_znfs.start):
            #         for z in allele_znfs.index:
            #             if end[i] >= allele_znfs.iloc[z].end:
            #                 read_znfs.append(allele_znfs.iloc[z].znfs)
            #             # if read ends alignment here, skip remaining znf iterations
            #             elif (end[i] > allele_znfs.iloc[z].start) & (end[i] < allele_znfs.iloc[z].end):
            #                 read_znfs.append(allele_znfs.iloc[z].znfs)
            #                 break

            # # reads that start alignment in znf region
            # elif start[i] >= flank_length:
            #     for z in allele_znfs.index:
            #         # skip over znfs that occur before the read starts aligning
            #         if start[i] > allele_znfs.iloc[z].end:
            #             continue
            #         if end[i] > allele_znfs.iloc[z].end:
            #             read_znfs.append(allele_znfs.iloc[z].znfs)
            #         # if read ends alignment here, skip remaining znf iterations
            #         elif (end[i] > allele_znfs.iloc[z].start) & (end[i] <= allele_znfs.iloc[z].end):
            #             read_znfs.append(allele_znfs.iloc[z].znfs)
            #             break
            #     # reads that end alignment in right flank region
            #     if end[i] > max(allele_znfs.end):
            #         read_flanks.append("right")
                    
            # if args.debug:
            #print("%s, %s, %d, %d, %s" % (read.name, best_allele, start, end, read_znfs), file=sys.stderr)
            # print("%s" % read.name, file=sys.stderr)
            # print(best_allele, file=sys.stderr)
            # print(start, file=sys.stderr)
            # print(end, file=sys.stderr)
            # print(read_znfs, file=sys.stderr)
            # print("---", file=sys.stderr)

        # add znf and flank counts, such that the sum of read_znf_counts + read_flank_counts = num_quality_reads
        # length of best_alleles already accounted for in length of read_znfs and read_flanks
        # read_flanks_znfs = read_znfs + read_flanks
        #print(read_znfs, file=sys.stderr)
        #print(read_flanks, file=sys.stderr)
        #print(read_flanks_znfs, file=sys.stderr)
        # for z in read_znfs:
        #     read_znf_counts[z] = read_znf_counts.get(z, 0) + (1/len(read_flanks_znfs))
        # for f in read_flanks:
        #     read_flank_counts[f] = read_flank_counts.get(f, 0) + (1/len(read_flanks_znfs))

        # if len(read_znfs) > 0:
        #     all_best_read_znfs.append(len(read_znfs)/len(best_allele))


        # if len(all_best_read_znfs) > 0:
        #     num_znfs_in_reads.append(sum(all_best_read_znfs)/len(all_best_read_znfs))
        # print("summary", file=sys.stderr)
        # print(sum(all_best_read_znfs), file=sys.stderr)
        # print(len(all_best_read_znfs), file=sys.stderr)
        # print(all_best_read_znfs, file=sys.stderr)
        # print(num_znfs_in_reads, file=sys.stderr)
        # print("---", file=sys.stderr)

# add 0 counts for znfs not observed in reads
# for i in znfs_list:
#     read_znf_counts[i] = read_znf_counts.get(i, 0)

# get mean number of znfs in reads to use as lambda
#mean_lambda = sum(num_znfs_in_reads)/len(num_znfs_in_reads)
#manual_lambda = args.lam

# additional debugging printing for znf counts
print("Number of quality reads: %d" % num_quality_reads, file=sys.stderr)
#print("Lambda based on mean znf counts in reads: %f" % mean_lambda, file=sys.stderr)
#print("Manually requested lambda: %f" % manual_lambda, file=sys.stderr)
# if args.debug:
#     print("znf and flank counts", file=sys.stderr)
#     for z, count in read_znf_counts.items():
#         print(z, '\t', count, file=sys.stderr)
#     for f, count in read_flank_counts.items():
#         print(f, '\t', count, file=sys.stderr)
#     #print(num_znfs_in_reads, file=sys.stderr)
print("znf and flank counts:", file=sys.stderr)
# for z, count in read_znf_counts.items():
#     print(z, '\t', count, file=sys.stderr)
for f, count in read_flank_counts.items():
    print(f, '\t', count, file=sys.stderr)


# get znf count profiles for each allele and compare to read count profiles:
# allele_sizes = znfs.groupby(['allele']).size().to_frame('allelesize').reset_index()
# all_allele_znf_counts = znfs[['allele', 'znfs']].groupby(['allele', 'znfs'])['znfs'].count().to_frame('count').reset_index()
# all_allele_znf_proportions = all_allele_znf_counts.merge(allele_sizes, how='outer')
# all_allele_znf_proportions['proportion'] = (all_allele_znf_proportions['count'] / all_allele_znf_proportions['allelesize'])

# for a in alleles:
#     allele_score = 0
#     allele_znf_counts = all_allele_znf_proportions[all_allele_znf_proportions.allele == a.name].set_index('znfs')
#     allele_znf_list = allele_znf_counts.index.tolist()

#     # calculate lambda
#     num_znfs = len(allele_znf_list)
#     calculated_lambda = calculate_lambda(num_quality_reads, num_znfs, read_length, znf_length)
#     print("Allele %s has a lambda of %f" % (a.name, calculated_lambda), file=sys.stderr)

#     for z in znfs_list:
#         if z not in allele_znf_list:
#             allele_count = 0
#         else:
#             allele_count = allele_znf_counts.loc[z, 'proportion']
#         if z not in read_znf_counts.keys():
#             read_count = 0
#         else:
#             read_count = read_znf_counts[z]
#         if allele_count == 0:
#             allele_score += log_poisson_pmf_manual(number=read_count, lam=0.005)
#             #print("%s, %s, %f, %f, %f, %f" % (a.name, z, read_count, allele_count, poisson.logpmf(read_count, 0.005), allele_score))
#         else:
#             # if manual_lambda > 0:
#             #     allele_score += log_poisson_pmf_manual(read_count, allele_count * manual_lambda)
#             #     #print("%s, %s, %f, %f, %f, %f" % (a.name, z, read_count, allele_count, poisson.logpmf(read_count, allele_count * manual_lambda), allele_score))
#             # else:
#             allele_score += log_poisson_pmf_manual(number=read_count, lam=(allele_count * calculated_lambda))
#                 #print("%s, %s, %f, %f, %f, %f" % (a.name, z, read_count, allele_count, poisson.logpmf(read_count, allele_count * mean_lambda), allele_score))
#     print("%s\t%f" % (a.name, allele_score))