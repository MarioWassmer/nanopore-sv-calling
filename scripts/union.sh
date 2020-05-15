#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: ./union.sh samples"
	exit 1
else

	mkdir -p union-set

	IN=$1
	OUT=union-set/union-set.vcf.gz

	if [[ "$IN" == *gz ]]
	then
		# debug
		echo "BGZIP compressed file detected."
	else
		# create union set using survivor
		SURVIVOR merge $IN 250 1 -1 -1 -1 49 union-set/union.vcf

		# sort, bgzip and tabix the union set
		cat <(cat union-set/union.vcf | grep "^#") <(cat union-set/union.vcf | grep -vE "^#" | sort -k1,1 -k2,2n) | bgzip -c > $OUT

		tabix $OUT

		rm union-set/union.vcf

		# call truvari on the union set
		truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/deletions-50-20k.vcf.gz -c $OUT --includebed /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed -r 1000 -p 0.00 -o truvari-union
	fi
fi
