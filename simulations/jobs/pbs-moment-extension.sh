#!/bin/sh
#PBS -N bfpwr-moment
#PBS -l select=1:ncpus=1:mem=2gb:scratch_local=1gb
#PBS -l walltime=06:00:00
#PBS -j oe
set -eu
: "${PHASE:?Set PHASE to plan, import, chunks, summaries, finalize, or assemble.}"
trap 'clean_scratch' TERM EXIT
SOURCE=${SOURCE:-/auto/brno2/home/fbartos/jobs/bfpwr/pr8-moment-extension}
INPUT=${INPUT:-/auto/brno2/home/fbartos/jobs/bfpwr/simulations/corpus/v1}
VERSION=${VERSION:-moment-v2}
OUTPUT=${OUTPUT:-/auto/brno2/home/fbartos/jobs/bfpwr/simulations/corpus/$VERSION}
RELEASE=${RELEASE:-"$SOURCE/release-$VERSION"}
MODULE_INIT=/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/libexec/modulecmd.tcl
MODULE_TCL=/cvmfs/software.metacentrum.cz/modulefiles/5.3.1/bin/tclsh
eval "$("$MODULE_TCL" "$MODULE_INIT" sh autoinit)"
module load r/4.5.1-gcc-10.2.1-zmneq6c
export R_LIBS_USER=${R_LIBS_USER:-/auto/brno2/home/fbartos/Rpackages-45}
export BFPWR_SIM_SOURCE_REVISION=${BFPWR_SIM_SOURCE_REVISION:-$(git -C "$SOURCE" rev-parse HEAD)}
cp -R "$SOURCE/package" "$SOURCE/simulations" "$SCRATCHDIR/"
cd "$SCRATCHDIR"
if [ "$PHASE" = import ]; then
    set -- simulations/scripts/run_moment_extension.R --phase import \
        --input-corpus "$INPUT" --output "$OUTPUT" --reuse-from "$REUSE_FROM" \
        --bf-source-root "$BF_SOURCE_ROOT"
    if [ -n "${OFFICIAL_METADATA:-}" ]; then
        set -- "$@" --official-metadata "$OFFICIAL_METADATA"
    fi
    Rscript "$@"
    exit 0
fi
Rscript simulations/scripts/run_moment_extension.R --phase "$PHASE" \
    --input-corpus "$INPUT" --output "$OUTPUT" --workers "${WORKERS:-1}"
if [ "$PHASE" = assemble ]; then
    Rscript simulations/scripts/verify_moment_search_limitation.R --output "$OUTPUT"
    Rscript simulations/scripts/validate_moment_extension.R --output "$OUTPUT"
    if [ "${PACKAGE:-true}" = true ]; then
        Rscript simulations/scripts/package_moment_extension.R \
            --output "$OUTPUT" --release "$RELEASE" --version "$VERSION"
    fi
fi
