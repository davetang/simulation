#!/usr/bin/env python3
#
# written by Dave Tang
#

import sys
import os.path
import argparse
# https://pysam.readthedocs.io/en/latest/api.html
import pysam

parser = argparse.ArgumentParser(description =
    """
    Subtract one BAM file from another BAM file using a read's
    query name.""",
    formatter_class = argparse.RawTextHelpFormatter
)
parser.add_argument(
    "bam_out",
    help = "Output BAM file",
)
parser.add_argument(
    "bam_whole",
    help = "BAM file to subtract from",
)
parser.add_argument(
    "bam_subset",
    help = "BAM file used for subtracting",
)
parser.add_argument(
    "-v",
    "--verbose",
    help = "verbose mode (default = False)",
    default = False,
    action = "store_true"
)
parser.add_argument(
        "-p",
        "--threads",
        help = "number of threads (default = 4)",
        choices = range(1,9),
        default = 4,
        type = int
)
args = parser.parse_args()

if os.path.isfile(args.bam_out):
    print(f"Specified BAM output {args.bam_out} already exists", file = sys.stderr)
    what_to_do = input("Continue y/N? ")
    if what_to_do == "y" or what_to_do == "Y":
        print("Continuing", file = sys.stderr)
    else:
        sys.exit()
print("Using {} threads".format(args.threads), file = sys.stderr)

if args.verbose:
    print("Verbose mode", file = sys.stderr)
    print("BAM whole: {}".format(args.bam_whole), file = sys.stderr)
    print("BAM subset: {}".format(args.bam_subset), file = sys.stderr)

def count_reads(bam):
    b = pysam.AlignmentFile(bam, "rb", threads = args.threads)
    print(f"Counting number of reads in {bam}", file = sys.stderr)
    entries = b.count(until_eof = True)
    print(f"{bam} has {entries} reads", file = sys.stderr)
    b.close()
    return(entries)

def which_pair(flag):
    if flag & 0x40:
        return(0)
    elif read.flag & 0x80:
        return(1)
    else:
        print("Could not determine read pair", file = sys.stderr)
        sys.exit()

entries = count_reads(args.bam_subset)
one_pc = round(entries*0.01) - 1
reads_processed = 0

reads = [{}, {}]
bam_sub = pysam.AlignmentFile(args.bam_subset, "rb", threads = args.threads)
for read in bam_sub.fetch(until_eof = True):
    read_pair = which_pair(read.flag)
    reads[read_pair][read.query_name] = 1
    reads_processed += 1
    if reads_processed % one_pc == 0:
        print(f'\rStoring {args.bam_subset}: {round(reads_processed / entries * 100)}% complete', end='', file = sys.stderr)
print(file = sys.stderr)
bam_sub.close()

entries = count_reads(args.bam_whole)
one_pc = round(entries*0.01) - 1
reads_processed = 0

bam_whole = pysam.AlignmentFile(args.bam_whole, "rb", threads = args.threads)
bam_out = pysam.AlignmentFile(args.bam_out, "wb", template = bam_whole)
for read in bam_whole.fetch(until_eof = True):
    read_pair = which_pair(read.flag)
    reads_processed += 1
    if reads_processed % one_pc == 0:
        print(f'\rProcessing {args.bam_whole}: {round(reads_processed / entries * 100)}% complete', end='', file = sys.stderr)

    if read.query_name in reads[read_pair]:
        continue
    else:
        bam_out.write(read)
print(file = sys.stderr)
print("Done", file = sys.stderr)

bam_whole.close()
bam_out.close()
sys.exit()

