#!/bin/bash
#
cp ../scripts/*.java .
docker build --progress=plain -t fcunial/callset_integration .
docker push fcunial/callset_integration
rm -f *.java *.class