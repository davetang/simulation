#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=samtools-${ver}
check_dir ${d}

wget https://github.com/samtools/samtools/releases/download/${ver}/samtools-${ver}.tar.bz2
tar xjf samtools-${ver}.tar.bz2
rm samtools-${ver}.tar.bz2

cd samtools-${ver}
./configure --prefix=$(pwd)
make && make install
cd ..
ln -v -s -f ${d}/samtools

bye
