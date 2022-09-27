#!/usr/bin/env python3
#
# written by Dave Tang
#

import sys
import argparse
from collections import defaultdict
# https://pysam.readthedocs.io/en/latest/api.html
import pysam

parser = argparse.ArgumentParser(description =
    """
    Read count summary of mapped dwgsim generated reads grouped by
    read pair, variant type, and number of sequencing errors""",
    formatter_class = argparse.RawTextHelpFormatter
)
parser.add_argument(
    "bam",
    help = "mapping results in BAM format of dwgsim generated reads",
)
parser.add_argument(
    "-v",
    "--verbose",
    help = "verbose mode (default = False)",
    default = False,
    action = "store_true"
)
parser.add_argument(
    "-u",
    "--unmapped",
    help = "include unmapped reads (default = False)",
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

if args.verbose:
    print("Verbose mode", file = sys.stderr)
    print("Input: {}".format(args.bam), file = sys.stderr)

print("Using {} threads".format(args.threads), file = sys.stderr)

def count_reads(bam):
    b = pysam.AlignmentFile(bam, "rb", threads = args.threads)
    print(f"Counting number of reads in {bam}", file = sys.stderr)
    entries = b.count(until_eof = args.unmapped)
    print(f"{bam} has {entries} reads", file = sys.stderr)
    b.close()
    return(entries)

entries = count_reads(args.bam)
one_pc = round(entries*0.01) - 1
reads_processed = 0

bam = pysam.AlignmentFile(args.bam, "rb", threads = args.threads)

stats = [{}, {}, {}]
for read in bam.fetch(until_eof = args.unmapped):

    if read.flag & 0x40:
        read_pair = 0
    elif read.flag & 0x80:
        read_pair = 1
    else:
        print("Could not determine read pair", file = sys.stderr)
        sys.exit()

    my_split = read.query_name.replace(":", "_").split("_")

    # contigs can have underscores, which causes a variable number of splits
    assert len(my_split) >= 14, my_split
    # assign backwards
    my_split.reverse()
    my_split[1:13] = [ int(x) for x in my_split[1:13] ]
    read_id, indel2, snp2, seq_error2, indel1, snp1, seq_error1, random2, \
    random1, strand2, strand1, start2, start1 = my_split[0:13]
    contig = '_'.join(my_split[13:])

    if read_pair == 0 and indel1 > 0 and snp1 > 0:
        var = 'both'
    elif read_pair == 1 and indel2 > 0 and snp2 > 0:
        var = 'both'
    elif indel1 > 0 or indel2 > 0:
        var = 'indel'
    elif snp1 > 0 or snp2 > 0:
        var = 'snp'
    else:
        var = 'none'

    if read_pair == 0:
        seq_error = seq_error1
    else:
        seq_error = seq_error2

    if var in stats[read_pair]:
        if seq_error in stats[read_pair][var]:
            stats[read_pair][var][seq_error] += 1
        else:
            stats[read_pair][var][seq_error] = 1
    else:
        stats[read_pair][var] = defaultdict(dict)
        stats[read_pair][var][seq_error] = 1

    reads_processed += 1
    if reads_processed % one_pc == 0:
        print(f'\rProcessing {args.bam}: {round(reads_processed / entries * 100)}% complete', end='', file = sys.stderr)
print(file = sys.stderr)

print("\t".join(["read", "var", "seq_error", "value"]))
for idx, read_pair in enumerate(stats):
    for var, seq_error in sorted(read_pair.items()):
        for seq_error, value in sorted(seq_error.items()):
            print(f"read{idx+1}\t{var}\t{seq_error}\t{value}")

print("Done", file = sys.stderr)
bam.close()
sys.exit()

