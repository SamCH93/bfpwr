#!/bin/sh
#PBS -N bfpwr-moment-archive
#PBS -l select=1:ncpus=1:mem=2gb:scratch_local=1gb
#PBS -l walltime=01:00:00
#PBS -j oe
set -eu
trap 'clean_scratch' TERM EXIT
SOURCE=${SOURCE:-/auto/brno2/home/fbartos/jobs/bfpwr/pr8-moment-extension}
export SOURCE
export VERSION=${VERSION:-moment-v2}
export OUTPUT=${OUTPUT:-/auto/brno2/home/fbartos/jobs/bfpwr/simulations/corpus/$VERSION}
export RELEASE=${RELEASE:-"$SOURCE/release-$VERSION"}
MODULE_INIT=/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/libexec/modulecmd.tcl
MODULE_TCL=/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/bin/tclsh
eval "$("$MODULE_TCL" "$MODULE_INIT" sh autoinit)"
module load r/4.5.1-gcc-10.2.1-zmneq6c
export R_LIBS_USER=${R_LIBS_USER:-/auto/brno2/home/fbartos/Rpackages-45}
# The common packager references original BF archives when BF import is recorded.
sh "$SOURCE/simulations/jobs/package-moment-extension.sh"
