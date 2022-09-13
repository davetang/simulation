#!/usr/bin/env bash

set -euo pipefail

threads=8
>&2 echo Using ${threads} threads

bin=$(dirname $0)
root=${bin}/..

ref=hg38
refdir=${root}/data/${ref}_dragmap

read1=${root}/test/${ref}.bwa.read1.fastq.gz
read2=${root}/test/${ref}.bwa.read2.fastq.gz

if [[ ! -e ${read1} || ! -e ${read2} ]]; then
   >&2 echo please run ${bin}/sim.sh to simulate reads
   exit 1
fi

${bin}/dragen-os \
   --num-threads ${threads} \
   -r ${refdir} \
   -1 ${read1} \
   -2 ${read2} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${root}/test/${ref}.dr.bam |\
   ${bin}/samtools index - ${root}/test/${ref}.dr.bam.bai

>&2 echo Done
exit 0

