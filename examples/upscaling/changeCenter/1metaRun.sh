#!/bin/bash

script="./2run.sh"

action=(mesh periodic relaxation)

for ((j = 100 ; j <= 200 ; j=j+25));
do
	for ((i = 0 ; i <= 200 ; i=i+25));
	do
		for a in ${action[@]}
		do
			echo "$a $i $j"
			$script $i $j $a > out${a}_${i}_${j}.txt
		done
	done
done

