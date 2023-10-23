#!/bin/bash
#
docker build --progress=plain -t fcunial/callset_integration .
docker push fcunial/callset_integration
