#!/bin/bash

cd ${0%/*} || exit 1    # Run from this directory

echo "Cleaning..."

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

clean_precice() {
  rm -fv ./precice-*-iterations.log \
    ./precice-*-convergence.log \
    ./precice-*-events.json \
    ./precice-*-events-summary.log \
    ./precice-postProcessingInfo.log \
    ./precice-*-watchpoint-*.log \
    ./precice-*-watchintegral-*.log \
    ./core
}

# Participant 1: Fluid
Participant1="Fluid"
(
  cd ${Participant1} || exit
  cleanCase0
  clean_precice
)

# Participant 2: Solid
Participant2="Solid"
(
  cd ${Participant2} || exit
  rm -fv *.log *.pla *.vtk *.out *.sock *.mov *.ine log.*
  clean_precice
)

# Remove the log files
rm -fv *.log
rm -fv log.*

# Remove the preCICE address files
clean_precice
rm -rfv precice-run
rm -fv .*.address

echo "Cleaning complete!"
#------------------------------------------------------------------------------
