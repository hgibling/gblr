#! /usr/bin/env python

### import libraries
import argparse
import os
import pysam
import re
import sys

### define functions
# get positions of deletions indels
def get_deletion_positions(read, region, gap_tolerance=args.gap_tolerance):
    split_cigar = [x for x in re.findall('[0-9]*[A-Z=]', read.cigarstring) if not 'S' in x]
    reference_pos = read.reference_start + 1
    deletion_positions = []

    # parse if cigar has deletions
    if any(chunk.endswith('D') for chunk in split_cigar):
        for chunk in split_cigar:
            # advance reference_pos for each match and deletion (insertions ignored as they are not in reference)
            if chunk.endswith(('M', 'X', '=')):
                print("%s match at %d" % (read.query_name, reference_pos))
                reference_pos += int(chunk[:chunk.find('MX=')])
            elif chunk.endswith('D'):
                # only consider deletions > 10bp
                if int(chunk[:chunk.find('D')]) > gap_tolerance:
                    deletion_positions.append(":".join([str(reference_pos), str(reference_pos + int(chunk[:chunk.find('D')]) - 1)]))
                    print("big del!")
                print("%s del at %d" % (read.query_name, reference_pos))
                reference_pos += int(chunk[:chunk.find('D')])
        if len(deletion_positions) == 0:
            deletion_positions.append('None')
    else:
        deletion_positions = ['None']
    return(deletion_positions)

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', type=str, required=True, help='bam file of aligned sequencing reads')
parser.add_argument('-r', '--region', type=str, default="5:23526782,23527873", help='position of one region of interest (ex: chr1:100000-200000)')
parser.add_argument('-g', '--gap-tolerance', type=int, default=10, help='ignore gaps this size or smaller')
parser.add_argument('-v', '--verbose', action='store_true', help='print stats about reads to stderr')
parser.add_argument('-o', '--output-name', type=str, required=True, help='name of file to filtered bam')
args = parser.parse_args()

### check arguments
if args.gap_tolerance < 0:
    exit("ERROR: gap tolerance must be a positive integer. Check parameters.")
if len(re.split("[:,-]", args.region)) != 3:
    exit("ERROR: region must have a chromosome, start, and end position. Check parameters.")

### get read sequences
try:
    bam = pysam.AlignmentFile(args.bam, 'rb')
except:
    exit("ERROR: issue with the reads file. Is it a proper bam file?")

### define variables
all_deletion_positions_dict = {'None': []}
keep_reads = set()
discard_reads = set()
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
    keep_reads.add(read.query_name)

    ### filter reads tobam  consider only those that touch both flanks
    if read.reference_start < (region[1] - args.flank_tolerance) and read.reference_end > (region[2] + args.flank_tolerance):

        ### collect positions of deletions
        read_deletions = get_deletion_positions(read, region)
        for pos in read_deletions:
            if pos not in all_deletion_positions_dict:
                all_deletion_positions_dict[pos] = [read.query_name]
            else:
                all_deletion_positions_dict[pos].append(read.query_name)

### get set of reads that have a large deletion at a unique position
all_reads_length = len(keep_reads)
for gap, reads in all_deletion_positions_dict.items():
    if len(reads) == 1:
        discard_reads.add(reads[0])

for read in discard_reads:
    keep_reads.remove(read)

### parse bam name
bam_name = re.split("[.]", args.bam)[0]
keep_reads_name = "."join(bam_name, "keep-reads.txt")

### write read list to file
file = open(keep_reads_name, "W")
for read in keep_reads:
    file.write("%s\n" % read)
file.close()

# print read stats
if args.verbose:
    print("Total number of reads: %d" % (all_reads_length), file=sys.stderr)
    print("Number of reads kept: %d" % (len(keep_reads)), file=sys.stderr)