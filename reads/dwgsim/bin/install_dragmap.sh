#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=DRAGMAP-${ver}
check_dir ${d}

wget https://github.com/Illumina/DRAGMAP/archive/refs/tags/${ver}.tar.gz
tar xzf ${ver}.tar.gz
rm ${ver}.tar.gz

cd DRAGMAP-${ver}
# https://github.com/Illumina/DRAGMAP/issues/22
HAS_GTEST=0 make CFLAGS:=
cd ..
ln -v -s -f ${d}/build/release/dragen-os

bye
