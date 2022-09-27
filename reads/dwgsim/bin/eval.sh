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

bam=$1
base=$(basename ${bam} .bam)
outdir=$(dirname ${bam})

${bin}/dwgsim_eval ${bam} |\
   gzip > ${outdir}/${base}.eval.txt.gz

${bin}/dwgsim_eval -i ${bam} |\
   gzip > ${outdir}/${base}.indels.eval.txt.gz

${bin}/dwgsim_eval -p ${bam} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${outdir}/${base}.incorrect.bam |\
   ${bin}/samtools index - ${outdir}/${base}.incorrect.bam.bai

${bin}/dwgsim_eval -p -i ${bam} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${outdir}/${base}.incorrect.indels.bam |\
   ${bin}/samtools index - ${outdir}/${base}.incorrect.indels.bam.bai

${bin}/bam_subtract.py ${outdir}/${base}.correct.bam ${bam} ${outdir}/${base}.incorrect.bam
${bin}/samtools index ${outdir}/${base}.correct.bam

${bin}/bam_summary.py ${bam} > ${outdir}/${base}.tsv
${bin}/bam_summary.py ${outdir}/${base}.correct.bam > ${outdir}/${base}.correct.tsv
${bin}/bam_summary.py ${outdir}/${base}.incorrect.bam > ${outdir}/${base}.incorrect.tsv

>&2 echo Done
exit 0

