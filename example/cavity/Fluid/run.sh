#!/bin/bash

cd ${0%/*} || exit 1
. $WM_PROJECT_DIR/bin/tools/RunFunctions

restore0Dir
touch res.foam

blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1

../../tools/run-openfoam.sh "$@"