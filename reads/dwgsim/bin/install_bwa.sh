#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=bwa-${ver}
check_dir ${d}

wget https://github.com/lh3/bwa/archive/refs/tags/v${ver}.tar.gz
tar xzf v${ver}.tar.gz
rm v${ver}.tar.gz

cd bwa-${ver}
make
cd ..
ln -v -s -f ${d}/bwa

bye
