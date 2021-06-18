#!/bin/bash

set -e
. $WM_PROJECT_DIR/bin/tools/RunFunctions

parallel=false
verbose=false
threads=false
case_name="Fluid"

run-openfoam() {
  solver=$(getApplication)
  if [ $parallel = true ]; then
    decomposePar -force > ../log."$case_name".decomposePar 2>&1
    nproc=$(getNumberOfProcessors)
    if [ $threads = true ]; then
      mpirun --use-hwthread-cpus -np "$nproc" "$solver" -parallel > ../log."$case_name" 2>&1 &
    else
      mpirun -np "$nproc" "$solver" -paralllel > ../log."$case_name" 2>&1 &
    fi
    run_mode="parallel"
  else
    "$solver" > ../log."$case_name" 2>&1 &
    run_mode="serial"
  fi
  solver_pid=$!
  trap '[ -z $solver_pid ] || kill $solver_pid; exit 1' SIGHUP SIGINT SIGQUIT SIGTERM
  echo "Solver '$solver' started in $run_mode with PID: $solver_pid"
  wait -f $solver_pid
  if [ $parallel = true ]; then
    reconstructPar
  fi
}

usage() {
  echo "Usage: run-openfoam [ -p | --parallel ] [ -t | --threads ]
                    [ -v | --verbose ] [ -c | --case CASE ]"
  exit 2
}

opts=$(getopt -a -n run-openfoam -o ptvc: --long parallel,threads,verbose,case: -- "$@")
if [ "$?" != "0" ]; then usage; fi

eval set -- "$opts"
while :
do
  case "$1" in
    -p | --parallel)   parallel=true; shift ;;
    -t | --threads)    threads=true; shift ;;
    -v | --verbose)    verbose=true; shift ;;
    -c | --case)       case_name="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."; usage ;;
  esac
done

if [ $verbose = true ]; then
  echo "verbose     : $verbose"
  echo "parallel    : $parallel"
  echo "threads     : $threads"
  echo "case        : $case_name"
  set -x
fi

run-openfoam
