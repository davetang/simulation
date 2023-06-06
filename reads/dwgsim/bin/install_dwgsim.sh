#!/usr/bin/env bash
#
# https://github.com/nh13/DWGSIM/blob/main/docs/02_Installation.md
#

set -euo pipefail

script_dir=$(dirname $0)
source ${script_dir}/utils.sh

check_args $# 1
ver=$1

cd ${script_dir}
d=dwgsim-${ver}
check_dir ${d}

git clone --recursive https://github.com/nh13/DWGSIM.git dwgsim-${ver}
cd dwgsim-${ver}
make clean
make
make test
cd ..
ln -v -s -f ${d}/dwgsim
ln -v -s -f ${d}/dwgsim_eval

bye
