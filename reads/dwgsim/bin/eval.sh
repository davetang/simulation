#!/usr/bin/env bash

set -euo pipefail

usage(){
   >&2 echo "Usage: $0 <infile.bam>"
   exit 1
}

if [[ $# < 1 ]]; then
   usage
fi

threads=8
bin=$(dirname $0)
root=${bin}/..

bam=$1
base=$(basename ${bam} .bam)

${bin}/dwgsim_eval ${bam} |\
   gzip > ${root}/test/${base}.eval.txt.gz

${bin}/dwgsim_eval -i ${bam} |\
   gzip > ${root}/test/${base}.indels.eval.txt.gz

${bin}/dwgsim_eval -p ${bam} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${root}/test/${base}.incorrect.bam |\
   ${bin}/samtools index - ${root}/test/${base}.incorrect.bam.bai

${bin}/dwgsim_eval -p -i ${bam} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${root}/test/${base}.incorrect.indels.bam |\
   ${bin}/samtools index - ${root}/test/${base}.incorrect.indels.bam.bai

>&2 echo Done
exit 0

