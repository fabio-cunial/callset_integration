#!/bin/bash
#
TAG=""
cp ../scripts/*.java .
docker build --progress=plain -t fcunial/callset_integration .
docker push fcunial/callset_integration${TAG}
rm -f *.java *.class