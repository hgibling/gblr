#! /usr/bin/env python

### import libraries
import argparse
import os
import pysam
import re
import sys
import warnings

### define functions
# get positions of indels
def get_indel_positions(read, region, gap_tolerance):
    split_cigar = [x for x in re.findall('[0-9]*[A-Z=]', read.cigarstring) if not 'S' in x]
    reference_pos = read.reference_start + 1
    deletion_positions = []
    insertion_positions = []

    # parse if cigar has indels
    if any(chunk.endswith(('D', 'I')) for chunk in split_cigar):
        for chunk in split_cigar:
            # advance reference_pos for each match and deletion (insertions ignored as they are not in reference)
            if chunk.endswith(('M', 'X', '=')):
                reference_pos += int(chunk[:chunk.find('MX=')])
            elif chunk.endswith('D'):
                # only consider deletions > gap_tolerance
                if int(chunk[:chunk.find('D')]) > gap_tolerance:
                    deletion_positions.append(":".join([str(reference_pos), str(reference_pos + int(chunk[:chunk.find('D')]) - 1)]))
                reference_pos += int(chunk[:chunk.find('D')])
            elif chunk.endswith('I'):
                # only consider insertions > gap_tolerance
                if int(chunk[:chunk.find('I')]) > gap_tolerance:
                    insertion_positions.append(":".join([str(reference_pos), str(reference_pos + int(chunk[:chunk.find('I')]) - 1)]))
        if len(deletion_positions) == 0:
            deletion_positions.append('None')
        if len(insertion_positions) == 0:
            insertion_positions.append('None')
        
    else:
        deletion_positions = ['None']
        insertion_positions = ['None']
    return([deletion_positions, insertion_positions])

# toss indels too infrequent or incorrect size
def toss_indels(gap_positions_dict, read_number_threshold, length_multiple, length_tolerance):
    discard_reads_set = set()
    for gap, reads in gap_positions_dict.items():
        if 'None' in gap:
            continue
        elif len(reads) < read_number_threshold:
            # toss reads that have gap at infrequent positions
            discard_reads_set.update(reads)
        elif args.length_multiple > 0:
            # toss reads with gap not a multiple of expected gap length
            gap_split = gap.split(":")
            gap_size = int(gap_split[1]) - int(gap_split[0]) + 1
            gap_size_multiple_modulo = gap_size % length_multiple
            if gap_size_multiple_modulo == 0:
                continue
            elif (gap_size_multiple_modulo > length_multiple / 2) and (gap_size_multiple_modulo < length_multiple - length_tolerance):
                # gap is too small and is tossed
                discard_reads_set.update(reads)
            elif (gap_size_multiple_modulo < length_multiple / 2) and (gap_size_multiple_modulo > length_tolerance):
                # gap is too large and is tossed
                discard_reads_set.update(reads)      
    return(discard_reads_set) 


### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', type=str, required=True, help='bam file of aligned sequencing reads')
parser.add_argument('-r', '--region', type=str, default='chr5:23526673,23527764', help='position of one region of interest (ex: chr1:100000-200000)')
parser.add_argument('-f', '--flank-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the flanking regions outside the region of interest')
parser.add_argument('-g', '--gap-tolerance', type=int, default=10, help='ignore gaps this size or smaller')
parser.add_argument('-R', '--read-threshold', type=float, default=0.25, help='proportion of reads that must have a specific deletion in order to consider the deletion at that position ok')
parser.add_argument('-l', '--length-multiple', type=int, default=0, help='toss reads with gaps that do not correspond to expected repeat length X or multiple thereof (default: ignore flag)')
parser.add_argument('-L', '--length-tolerance', type=int, default=0, help='allow for length multiple to be within X bases (default: 0)')
parser.add_argument('-v', '--verbose', action='store_true', help='print stats about reads to stderr')
parser.add_argument('-o', '--output-name', type=str, default='', help='name of output file of reads to keep (default: BAMNAME-keep-reads.txt)')
args = parser.parse_args()

### check arguments
if args.gap_tolerance < 0:
    exit("ERROR: gap tolerance must be a positive integer. Check parameters.")
if args.read_threshold < 0 or args.read_threshold > 1:
    exit("ERROR: read threshold must be a proportion (between 0 and 1). Check parameters.")
if len(re.split("[:,-]", args.region)) != 3:
    exit("ERROR: region must have a chromosome, start, and end position. Check parameters.")
if args.length_multiple < 0:
    warnings.warn("Gap length filter should be a positive integer if intended for use. Negative integers and zero are ignored and no length filter will be applied.")
if args.length_tolerance < 0:
    warnings.warn("Gap length tolerance should be a positive integer if intended for use. Negative integers and zero are ignored and no length tolerance will be use (i.e., length multiples will only count if exact).")
if args.length_multiple < 1 & args.length_tolerance > 0:
    warnings.warn("Gap length-multiple filter tolerance not declared, so length tolerance value ignored.")

### get read sequences
try:
    bam = pysam.AlignmentFile(args.bam, 'rb')
except:
    exit("ERROR: issue with the reads file. Is it a proper bam file?")

### define variables
all_deletion_positions_dict = {'None': []}
all_insertion_positions_dict = {'None': []}
keep_reads = set()
region = re.split("[:,-]", args.region)
region = [region[0], int(region[1]), int(region[2])]

### check chromosome naming convention and update region if needed
if ("chr" in region[0]) != ("chr" in bam.references[0]):
    if "chr" in bam.references[0]:
        region[0] = "chr" + region[0]
    else:
        region[0] = region[0].replace("chr", "")

### iterate over reads that overlap with region of interest
for read in bam.fetch(region[0], region[1], region[2]):

    ### filter reads to consider only those that touch both flanks
    if read.reference_start < (region[1] - args.flank_tolerance) and read.reference_end > (region[2] + args.flank_tolerance):
        keep_reads.add(read.query_name)

        ### collect positions of indels
        read_indels = get_indel_positions(read, region, args.gap_tolerance)
        # deletions
        for pos in read_indels[0]:
            if pos not in all_deletion_positions_dict:
                all_deletion_positions_dict[pos] = [read.query_name]
            else:
                all_deletion_positions_dict[pos].append(read.query_name)
        # insertions
        for pos in read_indels[1]:
            if pos not in all_insertion_positions_dict:
                all_insertion_positions_dict[pos] = [read.query_name]
            else:
                all_insertion_positions_dict[pos].append(read.query_name)

### get set of reads to toss
all_reads_length = len(keep_reads)
read_number_threshold = round(all_reads_length * args.read_threshold)

discard_reads_deletions = toss_indels(all_deletion_positions_dict, read_number_threshold, args.length_multiple, args.length_tolerance)
discard_reads_insertions = toss_indels(all_insertion_positions_dict, read_number_threshold, args.length_multiple, args.length_tolerance)
discard_reads = discard_reads_deletions | discard_reads_insertions

### remove tossed reads from keep set
for read in discard_reads:
    keep_reads.remove(read)

### parse bam name
bam_name = re.sub(".bam", "", args.bam)
if args.output_name != "":
    keep_reads_name = args.output_name
else:
    keep_reads_name = ".".join([bam_name, "keep-reads.txt"])

### write keep read list to file
file = open(keep_reads_name, "w")
for read in keep_reads:
    file.write("%s\n" % read)
file.close()

# print read stats
if args.verbose:
    print("Total number of reads: %d" % (all_reads_length), file=sys.stderr)
    print("Number of reads kept: %d (%f%%)" % (len(keep_reads), round((len(keep_reads)/all_reads_length * 100), 2)), file=sys.stderr)