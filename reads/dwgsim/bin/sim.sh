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
# -e FLOAT per base/color/flow error rate of the first read [from 0.020 to 0.020 by 0.000]
# -E FLOAT per base/color/flow error rate of the second read [from 0.020 to 0.020 by 0.000]

# note that the mutation rate is with respect to the reference FASTA. The
# default rate of 0.0010 will generate around 3 million mutations in a human
# reference genome
# -r FLOAT rate of mutations [0.0010]

# -R FLOAT fraction of mutations that are indels [0.10]
# -X FLOAT probability an indel is extended [0.30]
# -I INT the minimum length indel [1]
# -y FLOAT probability of a random DNA read [0.05]
# -z INT random seed (-1 uses the current time) [-1]

${bin}/dwgsim \
   -N 10000000 \
   -1 100 \
   -2 100 \
   -e 0.05 \
   -E 0.05 \
   -R 0.20 \
   -y 0.05 \
   -z 1984 \
   ${fasta} \
   ${root}/test/${ref}

>&2 echo Done
exit 0

