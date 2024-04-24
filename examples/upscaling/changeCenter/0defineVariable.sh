#!/bin/bash

hn=$(hostname)

if [ "$hn" = "GSE35838" ]
then
	np=1
	clusterString=" "
fi

if [ "$hn" = "achilles.unil.ch" ]
then
	np=2
	npReserve=2
	node="-q highmem -m node15 "
	#node="-m node05 "
	jobName="$3_$1_$2"
	outName="$jobName.out"
	errName="$jobName.err"
	echo "jobName= "$jobName
	clusterString="bsub $node -n $npReserve,$npReserve -J ${jobName} -o ${outName} -e ${errName} "
fi

parrotString="mpirun -n ${np} ../../../parrotc-opt -i "

export parrotString
export clusterString

