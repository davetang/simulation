#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=bcftools-${ver}
check_dir ${d}

wget https://github.com/samtools/bcftools/releases/download/${ver}/bcftools-${ver}.tar.bz2
tar xjf bcftools-${ver}.tar.bz2
rm bcftools-${ver}.tar.bz2

cd bcftools-${ver}
./configure --prefix=$(pwd)
make && make install
cd ..
ln -v -s -f ${d}/bcftools

bye
