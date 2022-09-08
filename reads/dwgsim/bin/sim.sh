#!/usr/bin/env bash

set -euo pipefail

bin=$(dirname $0)
root=${bin}/..

fasta=${root}/data/hg38.fa
ref=$(basename ${fasta} .fa)

if [[ ! -e ${fasta} ]]; then
   >&2 echo please run ${root}/data/download.sh to prepare ${fasta}
   exit 1
fi

${bin}/dwgsim -N 10000 -1 100 -2 100 -y 0 ${fasta} ${root}/test/${ref}

>&2 echo Done
exit 0

