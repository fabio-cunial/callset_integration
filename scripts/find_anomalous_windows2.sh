#!/bin/bash
INPUT_DIR=$1

set -euo pipefail


find ${INPUT_DIR} -type d -mindepth 1 -maxdepth 1 -print > directories.txt
while read DIRECTORY; do
    java CoverageStats ${DIRECTORY}
done < directories.txt
