#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l SVMerger.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2Align.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2DownloadAssemblies.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2Download.wdl
