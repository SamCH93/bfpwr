#!/bin/sh
# Append new moment-v2 assets; retain unchanged v1 BF archive references.
set -eu
SOURCE=${SOURCE:-$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)}
OUTPUT=${OUTPUT:-"$SOURCE/simulations/corpus/moment-v2"}
RELEASE=${RELEASE:-"$SOURCE/simulations/downloads/sim-corpus-v1-moment-v2"}
VERSION=${VERSION:-moment-v2}
"${RSCRIPT:-Rscript}" "$SOURCE/simulations/scripts/package_moment_extension.R" \
    --output "$OUTPUT" --release "$RELEASE" --version "$VERSION"
