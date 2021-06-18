#!/bin/bash

cd ${0%/*} || exit 1

../../tools/run-mbdyn.sh "$@"
