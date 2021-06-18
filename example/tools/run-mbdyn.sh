#!/bin/bash

set -eu

verbose=false
case_name="Solid"
script_name="predyn.py"

run-mbdyn() {
  python3 "$script_name" > ../log."$case_name" 2>&1 &
  solver_pid=$!
  trap '[ -z $solver_pid ] || kill $solver_pid; exit 1' SIGHUP SIGINT SIGQUIT SIGTERM
  echo "MBDyn started by Script '$script_name' with PID: $solver_pid"
  wait -f $solver_pid
}

usage() {
  echo "Usage: run-mbdyn [ -v | --verbose ] [ -c | --case CASE ]
                 [ -s | --script  SCRIPT ]"
  exit 2
}

opts=$(getopt -a -n run-openfoam -o vc:s: --long verbose,case_name:,script_name: -- "$@")
if [ "$?" != "0" ]; then usage; fi

eval set -- "$opts"
while :
do
  case "$1" in
    -v | --verbose)    verbose=true; shift ;;
    -c | --case)       case_name="$2"; shift 2 ;;
    -s | --script)     script_name="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."; usage ;;
  esac
done

if [ $verbose = true ]; then
  echo "verbose     : $verbose"
  echo "script name : $script_name"
  echo "case        : $case_name"
  set -x
fi

run-mbdyn
