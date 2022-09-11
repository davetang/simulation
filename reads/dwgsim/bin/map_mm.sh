#!/usr/bin/env bash

set -euo pipefail

threads=8
>&2 echo Using ${threads} threads

bin=$(dirname $0)
root=${bin}/..

fasta=${root}/data/hg38.fa
ref=$(basename ${fasta} .fa)

read1=${root}/test/${ref}.bwa.read1.fastq.gz
read2=${root}/test/${ref}.bwa.read2.fastq.gz

if [[ ! -e ${read1} || ! -e ${read2} ]]; then
   >&2 echo please run ${bin}/sim.sh to simulate reads
   exit 1
fi

# short genomic paired-end reads
${bin}/minimap2 -t ${threads} -ax sr ${fasta} ${read1} ${read2} |\
   ${bin}/samtools sort -@ ${threads} -O BAM |\
   tee ${root}/test/${ref}.mm.bam |\
   ${bin}/samtools index - ${root}/test/${ref}.mm.bam.bai

>&2 echo Done
exit 0

