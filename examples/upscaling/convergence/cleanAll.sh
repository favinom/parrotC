#!/bin/bash

for i in `ls output*csv`;
do
	./cleanFile.sh $i
done
