#!/bin/bash

# number of background elements
xcBash=$1
# number of amr steps
ycBash=$2

source ./0defineVariable.sh

if [ "$3" = "mesh" ]
then
$clusterString $parrotString mooseScripts/meshgen.i xcParrot=${xcBash} ycParrot=${ycBash}
fi

if [ "$3" = "periodic" ]
then
$clusterString $parrotString mooseScripts/runPeriodic.i xcParrot=${xcBash} ycParrot=${ycBash}
fi

if [ "$3" = "relaxation" ]
then
$clusterString $parrotString mooseScripts/runRelaxation.i xcParrot=${xcBash} ycParrot=${ycBash}
fi
