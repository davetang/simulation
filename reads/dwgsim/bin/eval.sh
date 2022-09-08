#!/usr/bin/env bash

set -euo pipefail

bin=$(dirname $0)
root=${bin}/..

bam=${root}/test/hg38.bam
ref=$(basename ${bam} .bam)

if [[ ! -e ${bam} ]]; then
   >&2 echo please run ${root}/bin/map.sh to generate ${bam}
   exit 1
fi

${bin}/dwgsim_eval ${bam} |\
   gzip > ${root}/test/${ref}.eval.txt.gz

>&2 echo Done
exit 0

