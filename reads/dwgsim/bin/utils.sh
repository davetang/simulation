#!/usr/bin/env bash

set -euo pipefail

usage(){
   >&2 echo "Usage: $0 <version>"
   exit 1
}

check_args(){
   if [[ $1 -ne $2 ]]; then
      usage
   fi
}

check_dir(){
   d=$1
   if [[ -d ${d} ]]; then
      >&2 echo ${d} already exists and will be removed
      rm -rf ${d}
   fi
}

bye(){
   >&2 echo Done
   exit 0
}
