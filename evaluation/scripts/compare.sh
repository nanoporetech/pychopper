#!/bin/bash

D1=$1
N1=$2
D2=$3
N2=$4

OD="cmp_${N1}_${N2}"
rm -fr $OD; mkdir -p $OD

trim_range(){
	cut -d$"|" -f 2 $1 | grep -v Read | sort -k 1,1 > $2
}

J_READ_COV=$OD/read_cov.tsv
echo -e "Read\tReadCov:${N1}\tReadCov:${N2}" > $J_READ_COV
trim_range $D1/read_cov.tsv $D1/clean_read_cov.tsv
trim_range $D2/read_cov.tsv $D2/clean_read_cov.tsv
join -t $'\t' -j 1 "$D1/clean_read_cov.tsv" "$D2/clean_read_cov.tsv" >> $J_READ_COV
./scripts/regplot.py "ReadCov:${N1}" "ReadCov:${N2}" $J_READ_COV $OD/read_cov_reg.pdf


J_REF_COV=$OD/ref_cov.tsv
echo -e "Read\tRefCov:${N1}\tRefCov:${N2}" > $J_REF_COV
trim_range $D1/ref_cov.tsv $D1/clean_ref_cov.tsv
trim_range $D2/ref_cov.tsv $D2/clean_ref_cov.tsv
join -t $'\t' -j 1 "$D1/clean_ref_cov.tsv" "$D2/clean_ref_cov.tsv" >> $J_REF_COV
./scripts/regplot.py "RefCov:${N1}" "RefCov:${N2}" $J_REF_COV $OD/ref_cov_reg.pdf
