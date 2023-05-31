#!/usr/bin/env bash
#
# https://github.com/nh13/DWGSIM/blob/main/docs/02_Installation.md
#

set -euo pipefail

cd $(dirname $0)

d=DWGSIM

bye(){
   >&2 echo Done
   exit 0
}

if [[ -d ${d} ]]; then
   >&2 echo ${d} already exists
   if [[ -e ${d}/dwgsim && ${d}/dwgsim_eval ]]; then
      >&2 echo dwgsim and dwgsim_eval already compiled
      if [[ ! -L dwgsim ]]; then
         >&2 echo Creating symlink to dwgsim
         ln -v -s ${d}/dwgsim
      else
         >&2 echo dwgsim symlink already exists
      fi
      if [[ ! -L dwgsim_eval ]]; then
         >&2 echo Creating symlink to dwgsim_eval
         ln -v -s ${d}/dwgsim_eval
      else
         >&2 echo dwgsim_eval symlink already exists
      fi
      bye
   fi
else
   git clone --recursive https://github.com/nh13/DWGSIM.git
fi

rm -f dwgsim dwgsim_eval
cd DWGSIM
make clean
make
make test
cd ..
ln -v -s ${d}/dwgsim
ln -v -s ${d}/dwgsim_eval

bye
