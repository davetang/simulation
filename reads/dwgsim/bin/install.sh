#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

./install_minimap2.sh
./install_dwgsim.sh
./install_samtools.sh

