#! /usr/bin/env python

### import libraries
import argparse
import os
import pysam
import re
import sys
import warnings

### define functions
# get positions of deletions indels
def get_deletion_positions(read, region, gap_tolerance):
    split_cigar = [x for x in re.findall('[0-9]*[A-Z=]', read.cigarstring) if not 'S' in x]
    reference_pos = read.reference_start + 1
    deletion_positions = []

    # parse if cigar has deletions
    if any(chunk.endswith('D') for chunk in split_cigar):
        for chunk in split_cigar:
            # advance reference_pos for each match and deletion (insertions ignored as they are not in reference)
            if chunk.endswith(('M', 'X', '=')):
                reference_pos += int(chunk[:chunk.find('MX=')])
            elif chunk.endswith('D'):
                # only consider deletions > 10bp
                if int(chunk[:chunk.find('D')]) > gap_tolerance:
                    deletion_positions.append(":".join([str(reference_pos), str(reference_pos + int(chunk[:chunk.find('D')]) - 1)]))
                reference_pos += int(chunk[:chunk.find('D')])
        if len(deletion_positions) == 0:
            deletion_positions.append('None')
    else:
        deletion_positions = ['None']
    return(deletion_positions)

### parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', type=str, required=True, help='bam file of aligned sequencing reads')
parser.add_argument('-r', '--region', type=str, default='chr5:23526673,23527764', help='position of one region of interest (ex: chr1:100000-200000)')
parser.add_argument('-f', '--flank-tolerance', type=int, default=50, help='minimum number of bases to which a read must align in the flanking regions outside the region of interest')
parser.add_argument('-g', '--gap-tolerance', type=int, default=20, help='ignore gaps this size or smaller')
parser.add_argument('-R', '--read-threshold', type=float, default=0.25, help='proportion of reads that must have a specific deletion in order to consider the deletion at that position ok')
parser.add_argument('-l', '--length-multiple', type=int, default=-1, help='toss reads with gaps that do not correspond to expected repeat length X or multiple thereof (default: ignore flag)')
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
if args.length_multiple < -1 | args.length_multiple == 0:
    warnings.warn("Gap length filter should be a positive integer if intended for use. Negative integers and zero are ignored and no length filter will be applied.")

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

    ### filter reads to consider only those that touch both flanks
    if read.reference_start < (region[1] - args.flank_tolerance) and read.reference_end > (region[2] + args.flank_tolerance):
        keep_reads.add(read.query_name)

        ### collect positions of deletions
        read_deletions = get_deletion_positions(read, region, args.gap_tolerance)
        for pos in read_deletions:
            if pos not in all_deletion_positions_dict:
                all_deletion_positions_dict[pos] = [read.query_name]
            else:
                all_deletion_positions_dict[pos].append(read.query_name)

### get set of reads that have a large deletion at a unique position
all_reads_length = len(keep_reads)
read_number_threshold = round(all_reads_length * args.read_threshold)
print("dict is %i, num reads threshold is %i" % (len(all_deletion_positions_dict), read_number_threshold), file=sys.stderr)
for gap, reads in all_deletion_positions_dict.items():
    # print("gap is %s" % gap, file=sys.stderr)
    if 'None' in gap:
        continue
    elif len(reads) < read_number_threshold:
        discard_reads.update(reads)
        # print("gap is %s, len reads is %i" % (gap, len(reads)), file=sys.stderr)
    elif args.length_multiple > 0:
        #print("gap is %s" % gap, file=sys.stderr)
        # toss reads if gap is not expected repeat length or multiple thereof
        gap_split = gap.split(":")
        if (int(gap_split[1]) - int(gap_split[0])) % args.length_multiple != 0:
            print("bad gap is %i" % (int(gap_split[1]) - int(gap_split[0])), file=sys.stderr)
            discard_reads.update(reads)
        else:
            print("good gap is %i" % (int(gap_split[1]) - int(gap_split[0])), file=sys.stderr)

for read in discard_reads:
    keep_reads.remove(read)

### parse bam name
bam_name = re.sub(".bam", "", args.bam)
if args.output_name != "":
    keep_reads_name = args.output_name
else:
    keep_reads_name = ".".join([bam_name, "keep-reads.txt"])

### write read list to file
file = open(keep_reads_name, "w")
for read in keep_reads:
    file.write("%s\n" % read)
file.close()

# print read stats
if args.verbose:
    print("Total number of reads: %d" % (all_reads_length), file=sys.stderr)
    print("Number of reads kept: %d (%f%%)" % (len(keep_reads), round((len(keep_reads)/all_reads_length * 100), 2)), file=sys.stderr)