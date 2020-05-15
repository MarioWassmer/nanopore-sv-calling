#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: $0 samples"
	exit 1
else

	mkdir -p intersection-set

	IN=$1
	OUT=intersection-set/intersection-set.vcf.gz

	if [[ "$IN" == *gz ]]
	then
		# debug
		echo "BGZIP compressed file detected."
	else
		# create union set using survivor
		SURVIVOR merge $IN 250 2 -1 -1 -1 49 intersection-set/intersection.vcf

		# sort, bgzip and tabix the union set
		cat <(cat intersection-set/intersection.vcf | grep "^#") <(cat intersection-set/intersection.vcf | grep -vE "^#" | sort -k1,1 -k2,2n) | bgzip -c > $OUT

		tabix $OUT

		rm intersection-set/intersection.vcf

		# call truvari on the union set
		truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/deletions.vcf.gz \
		 -c $OUT --includebed \
		 /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
		 --passonly \
		 -s 50 \
		 -r 1000 \
		 -p 0.00 \
		 -o deletions-intersect
	fi
fi
