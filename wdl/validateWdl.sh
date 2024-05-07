#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeDipcallSimple.wdl
java -jar ${WOMTOOL_PATH} validate -l RepeatAnnotation.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfCutesv.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfLrcaller.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfSvjedigraph.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfKanpig.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfSniffles.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeIntrasample.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2Align2.wdl
java -jar ${WOMTOOL_PATH} validate -l Dipcall2SVs.wdl
java -jar ${WOMTOOL_PATH} validate -l Hapdiff.wdl
java -jar ${WOMTOOL_PATH} validate -l AnnotateJointVcf.wdl
java -jar ${WOMTOOL_PATH} validate -l AddTruvariAnnotations.wdl
java -jar ${WOMTOOL_PATH} validate -l MergeRegenotypedIntersampleVcf.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcfKanpigMerged.wdl
java -jar ${WOMTOOL_PATH} validate -l ToBed.wdl
java -jar ${WOMTOOL_PATH} validate -l ConcatenateChromosomes.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l FilterAndSplit.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntrasampleFabio.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedVcf.wdl
java -jar ${WOMTOOL_PATH} validate -l GetRegenotypedTruePositives.wdl
java -jar ${WOMTOOL_PATH} validate -l ROC.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntersampleNaive.wdl
java -jar ${WOMTOOL_PATH} validate -l GraphEvaluationFabio.wdl
java -jar ${WOMTOOL_PATH} validate -l JasmineIntersample2.wdl
java -jar ${WOMTOOL_PATH} validate -l SvpopIntra.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeDipcall.wdl
java -jar ${WOMTOOL_PATH} validate -l PbsvIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMerge.wdl
java -jar ${WOMTOOL_PATH} validate -l HPRCDownloadClones.wdl
java -jar ${WOMTOOL_PATH} validate -l SVMergerIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l SVMerger.wdl
java -jar ${WOMTOOL_PATH} validate -l SnifflesIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l JasmineIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l Jasmine.wdl
java -jar ${WOMTOOL_PATH} validate -l SvimmerIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l Svimmer.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2Align.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2DownloadAssemblies.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC2Download.wdl
