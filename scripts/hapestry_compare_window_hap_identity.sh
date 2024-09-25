#!/bin/bash
#
ROOT_DIR=$1

set -euo pipefail

ls ${ROOT_DIR} > list.txt
while read DIRECTORY; do
    if [ -d ${ROOT_DIR}/${DIRECTORY} ]; then
        java CoverageStats7 ${ROOT_DIR}/${DIRECTORY}
    else 
        echo "${ROOT_DIR}/${DIRECTORY} is not a directory"
    fi
done < list.txt
rm -f list.txt
