IN=$1
cut -d, -f1,26 ${IN} | sort| uniq| grep -v asv_na| grep -v BC10 > ${IN%%.csv}_asv_indels.csv
python make_matrix.py ${IN%%.csv}_asv_indels.csv > ${IN%%.csv}_matrix.csv
