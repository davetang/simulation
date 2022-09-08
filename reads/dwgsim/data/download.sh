#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

urls=(
   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
   https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
)

for url in ${urls[@]}; do

   genome=$(basename ${url})
   
   if [[ ! -e ${genome} ]]; then
      >&2 echo Download ${genome}
      wget -c ${url}
   else
      >&2 echo ${genome} already exists
   fi

   genome_unzipped=$(basename ${url} .gz)
   if [[ ! -e ${genome_unzipped} ]]; then
      >&2 echo Extracting ${genome}
      gunzip -c ${genome} > ${genome_unzipped}
   else
      >&2 echo ${genome_unzipped} already exists
   fi
   
done

>&2 echo Done
exit 0

