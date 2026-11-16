#!/bin/bash
# Launcher for the SVengers Nextflow pipeline on maestro/helix.
#   ./run_pipeline.sh -params-file params.yml [extra nextflow args]
# Add -resume to continue an interrupted run.

set -euo pipefail
cd "$(dirname "$0")"

module load nextflow 2>/dev/null || true
module load java 2>/dev/null || true

# non-interactive shells: fall back to explicit paths
if ! command -v nextflow >/dev/null; then
    export PATH=/opt/gensoft/exe/nextflow/26.04.2/bin:$PATH
fi
if ! java -version 2>/dev/null; then
    export JAVA_HOME=/opt/gensoft/exe/graalvm/ce-java21-21.0.2
    export PATH=$JAVA_HOME/bin:$PATH
fi

export NXF_OPTS="-Xms500m -Xmx4g"

nextflow run main.nf \
    -profile slurm \
    "$@"
