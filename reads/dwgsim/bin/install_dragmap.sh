#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

ver=1.3.0
tool=DRAGMAP
binary=dragen-os
d=${tool}-${ver}

bye(){
   >&2 echo Done
   exit 0
}

if [[ -d ${d} ]]; then
   >&2 echo ${d} already exists
   if [[ -e ${d}/build/release/${binary} ]]; then
      >&2 echo ${binary} already compiled
      if [[ ! -L ${binary} ]]; then
         >&2 echo Creating symlink to ${binary}
         ln -v -s ${d}/build/release/${binary} .
      else
         >&2 echo ${binary} symlink already exists
      fi
      bye
   fi
else
   wget https://github.com/Illumina/DRAGMAP/archive/refs/tags/${ver}.tar.gz
   tar xzf ${ver}.tar.gz
   rm ${ver}.tar.gz
fi

rm -f ${binary}
cd ${tool}-${ver}
HAS_GTEST=0 make
cd ..
ln -v -s ${d}/build/release/${binary} .

bye

