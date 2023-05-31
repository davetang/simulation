#!/usr/bin/env bash

docker run --rm -it -u $(stat -c "%u:%g" $HOME) -v $(pwd):$(pwd) -w $(pwd) davetang/build:1.2.6 /bin/bash
