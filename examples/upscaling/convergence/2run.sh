#!/bin/bash

# number of background elements
nelBash=$1
# number of amr steps
amrBash=$2
# do boundary refinement boundary r
dbrStringBash=$3

if [ "$dbrStringBash" = "f" ]
then
	echo 'doing full refinement'
	dbrBash=false
else
	if [ "$dbrStringBash" = "b" ]
	then
		echo 'doing boundary refinement'
		dbrBash=true
	else
		echo 'error'
		exit 1
	fi
fi

source ./0defineVariable.sh

if [ "$4" = "mesh" ]
then
$clusterString $parrotString mooseScripts/meshgen.i nel=${nelBash} amr=${amrBash} dbr=${dbrBash} dbrString=${dbrStringBash}
fi

if [ "$4" = "run" ]
then
$clusterString $parrotString mooseScripts/run.i nel=${nelBash} amr=${amrBash} dbrString=${dbrStringBash}
fi
