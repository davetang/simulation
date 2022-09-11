#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

ver=0.7.17
tool=bwa
d=${tool}-${ver}

bye(){
   >&2 echo Done
   exit 0
}

if [[ -d ${d} ]]; then
   >&2 echo ${d} already exists
   if [[ -e ${d}/${tool} ]]; then
      >&2 echo ${tool} already compiled
      if [[ ! -L ${tool} ]]; then
         >&2 echo Creating symlink to ${tool}
         ln -v -s ${d}/${tool}
      else
         >&2 echo ${tool} symlink already exists
      fi
      bye
   fi
else
   wget https://github.com/lh3/bwa/archive/refs/tags/v${ver}.tar.gz
   tar xzf v${ver}.tar.gz
   rm v${ver}.tar.gz
fi

rm -f ${tool}
cd ${tool}-${ver}
make
cd ..
ln -v -s ${d}/${tool}

bye

