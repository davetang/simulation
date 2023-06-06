#!/usr/bin/env bash
#
# https://github.com/lh3/minimap2
#

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=minimap2-${ver}
check_dir ${d}

wget https://github.com/lh3/minimap2/archive/refs/tags/v${ver}.tar.gz
tar xzf v${ver}.tar.gz
rm v${ver}.tar.gz

cd minimap2-${ver}
make clean
make
cd ..
ln -v -s -f ${d}/minimap2

bye
