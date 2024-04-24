#!/bin/bash

for i in {0..20..5}
do
  echo "Number: $i"
done

exit 1

script="./2run.sh"

action=(mesh periodic relaxation)
boundaryType=(b f)

for j in {0..200..25}
do
	echo "$j"
	for i in {0..200..25}
	do
		for a in ${action[@]}
		do
			
			#$script $i $j $a
		done
	done
done
