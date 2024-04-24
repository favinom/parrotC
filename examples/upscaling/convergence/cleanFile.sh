#!/bin/bash

filename="new$1"

cp $1 $filename

sed -i -e 's/(//g' $filename
sed -i -e 's/,-0)//g' $filename
sed -i -e 's/,0)//g' $filename
tail +3 $filename > temp.csv
cat temp.csv > $filename
touch temp.csv
rm temp.csv

