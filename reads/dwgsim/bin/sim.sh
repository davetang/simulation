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

# -N INT number of read pairs (-1 to disable) [-1]
# -1 INT length of the first read [70]
# -2 INT length of the second read [70]
# -y FLOAT probability of a random DNA read [0.05]
# -z INT random seed (-1 uses the current time) [-1]

${bin}/dwgsim \
   -N 10000 \
   -1 100 \
   -2 100 \
   -y 0 \
   -z 1984 \
   ${fasta} \
   ${root}/test/${ref}

>&2 echo Done
exit 0

