#!/usr/bin/env bash
#
# https://github.com/lh3/minimap2
#

set -euo pipefail

cd $(dirname $0)

minimap2_ver=2.26
d=minimap2-${minimap2_ver}

bye(){
   >&2 echo Done
   exit 0
}

if [[ -d ${d} ]]; then
   >&2 echo ${d} already exists
   if [[ -e ${d}/minimap2 ]]; then
      >&2 echo minimap2 already compiled
      if [[ ! -L minimap2 ]]; then
         >&2 echo Creating symlink to minimap2
         ln -v -s ${d}/dwgsim
      else
         >&2 echo minimap2 symlink already exists
      fi
      bye
   fi
else
   wget https://github.com/lh3/minimap2/archive/refs/tags/v${minimap2_ver}.tar.gz
   tar xzf v${minimap2_ver}.tar.gz
   rm v${minimap2_ver}.tar.gz
fi

rm -f minimap2
cd minimap2-${minimap2_ver}
make clean
make
cd ..
ln -v -s ${d}/minimap2

bye
