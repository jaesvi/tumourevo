#!/usr/bin/env bash
set -euo pipefail

cd /home/valleinclan/dev/pipelines/tumourevo

nextflow run . \
    -profile singularity \
    -params-file test-purple-pass01/params.yml \
    -work-dir test-purple-pass01/work \
    -log test-purple-pass01/.nextflow.log \
    -c test-purple-pass01/macondo.config \
    -resume
