#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

ver=1.17
d=samtools-${ver}

bye(){
   >&2 echo Done
   exit 0
}

if [[ -d ${d} ]]; then
   >&2 echo ${d} already exists
   if [[ -e ${d}/samtools ]]; then
      >&2 echo samtools already compiled
      if [[ ! -L samtools ]]; then
         >&2 echo Creating symlink to samtools
         ln -v -s ${d}/samtools
      else
         >&2 echo samtools symlink already exists
      fi
      bye
   fi
else
   wget https://github.com/samtools/samtools/releases/download/${ver}/samtools-${ver}.tar.bz2
   tar xjf samtools-${ver}.tar.bz2
   rm samtools-${ver}.tar.bz2
fi

rm -f samtools
cd samtools-${ver}
./configure --prefix=$(pwd)
make && make install
cd ..
ln -v -s ${d}/samtools

bye
