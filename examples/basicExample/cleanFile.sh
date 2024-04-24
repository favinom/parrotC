sed -i -e 's/(//g' $1
sed -i -e 's/,-0)//g' $1
sed -i -e 's/,0)//g' $1
tail +3 $1 > temp.csv
cat temp.csv > new$1
touch temp.csv
rm temp.csv

