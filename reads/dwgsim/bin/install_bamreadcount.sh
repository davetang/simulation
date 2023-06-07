#!/usr/bin/env bash

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=bam-readcount-${ver}
check_dir ${d}

wget https://github.com/genome/bam-readcount/archive/v${ver}.tar.gz -O ${d}.tar.gz
tar xzf ${d}.tar.gz
rm ${d}.tar.gz

cd ${d}
mkdir build && cd build
cmake ..
export CFLAGS="-fPIC"
make

cd ../..
ln -v -s -f ${d}/build/bin/bam-readcount

bye
