#!/usr/bin/env bash

set -euo pipefail

if [[ ! -x $(command -v docker) ]]; then
   >&2 echo Could not find docker
   exit 1
fi

cd $(dirname $0)/..

docker run \
   --rm \
   -it \
   -u $(id -u):$(id -g) \
   -v $(pwd):/work \
   davetang/build:1.2.2 \
   /bin/bash

