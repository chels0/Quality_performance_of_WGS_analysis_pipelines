csplit fastqc_data.txt '/>>/' '{*}'
awk '{ print $1,$2 }' xx03 > xx03_relevant

tail -n +2 xx03_relevant > xx03.csv
