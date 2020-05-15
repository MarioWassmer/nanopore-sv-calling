#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage: ./survivor-eval.sh input.vcf"
	exit 1
else
	IN=$1
	FOLDER="$(dirname $(readlink -f IN))/SurvivorStats"
	mkdir -p $FOLDER
	if [[ "$1" != *gz ]]
	then
		printf "BGZIP input file...."
		bgzip $IN
		printf "DONE!\n"
	fi
	printf "Prepare input file...."
	cat <(zcat $1 | grep "^#") \
	<(zcat $1 | grep -vE "^#" | grep -v "GL000" \
	| grep -v "hs37" | grep -v "NC" | grep -v "MT" \
	| grep -w "PASS") > ${FOLDER}/survivor-prepared.vcf
	printf "DONE!\n"

	printf "Compute statistics...."
	SURVIVOR stats ${FOLDER}/survivor-prepared.vcf 50 20000 -1 ${FOLDER}/survivor-stats.txt > ${FOLDER}/survivor-summary.txt
	printf "DONE!\n\n"
	cat ${FOLDER}/survivor-summary.txt
fi
