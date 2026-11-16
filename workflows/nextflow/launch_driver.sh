#!/bin/bash
# Detached Nextflow driver for SVengers cohort runs on maestro/helix.
#
# The driver must survive the launching session (a driver killed mid-run
# leaves unrecorded task completions and stale lock files), so launch it
# detached:
#
#   setsid nohup ./launch_driver.sh \
#       --manifest .../COHORT.manifest.tsv \
#       --pedigree .../COHORT.pedigree.tsv \
#       --outdir .../runs/COHORT \
#       --cohort_name COHORT \
#       --steps call,filter,merge,genotype,cohort,ped,report \
#       -with-report -with-trace -resume \
#       > driver_COHORT.log 2>&1 < /dev/null & disown
#
# Forces the graalvm Java (the `java` module can resolve to an unsupported
# version in non-interactive shells) and pins the Nextflow version.
set -euo pipefail
cd "$(dirname "$0")"

export JAVA_HOME=/opt/gensoft/exe/graalvm/ce-java21-21.0.2
export JAVA_CMD=$JAVA_HOME/bin/java
export PATH=$JAVA_HOME/bin:/opt/gensoft/exe/nextflow/26.04.2/bin:$PATH
export NXF_OPTS="-Xms500m -Xmx4g"

exec nextflow run main.nf -profile slurm "$@"
