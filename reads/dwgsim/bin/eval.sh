#!/usr/bin/env bash

set -euo pipefail

usage(){
   >&2 echo "Usage: $0 <infile.bam>"
   exit 1
}

if [[ $# < 1 ]]; then
   usage
fi

bin=$(dirname $0)
root=${bin}/..

bam=$1
base=$(basename ${bam} .bam)

${bin}/dwgsim_eval ${bam} |\
   gzip > ${root}/test/${base}.eval.txt.gz

${bin}/dwgsim_eval -i ${bam} |\
   gzip > ${root}/test/${base}.indels.eval.txt.gz

>&2 echo Done
exit 0

