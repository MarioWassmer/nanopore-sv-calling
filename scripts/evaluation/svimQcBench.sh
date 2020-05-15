#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "Usage: $0 sniffles_output.vcf[.gz] threshold(int in [1,100])"
	exit 1
else
	# Save input file in variable
	IN=$1
	THRESHOLD=$2

	# get absolute path of this script and python scripts
	SCRIPTPATH="$(dirname $(readlink -f $0))"

	# get absolute path of folder where the script has been started
	FOLDER="$(dirname $(readlink -f IN))/qcBenchmark"
	mkdir -p $FOLDER

	# define truvari prepared file
	OUT=$FOLDER/$(echo $(basename $IN) | cut -f 1 -d '.')_ins-del.vcf.gz

	# define deletions and insertions output
	DEL=$FOLDER/deletions.vcf.gz
	INS=$FOLDER/insertions.vcf.gz

	# truvari output folder
	TRUV_DEL=$FOLDER/truvari-deletions
	TRUV_INS=$FOLDER/truvari-insertions

	# check if input file exists
	if [ ! -f "$IN" ]; then
		echo "$IN not found!"
		exit 1
	fi



	# bgzip file if it is not already
	if [[ "$IN" != *gz ]]
	then
		printf "Compress and index input file...."

		bgzip $IN

		IN=$(echo $(basename $IN).gz)

		printf "DONE!\n"
	fi

	# ================ Prepare for Truvari =================
	printf "Prepare input file for truvari evaluation...."
	# fix the header of the input file
	# Extract Insertions, Deletions and Duplications; Normalize Duplications to Insertions
	cat <(zcat $IN | grep "^#") \
        <(zcat $IN | grep -vE "^#" \
	| egrep "^[(0-9)(X)(Y)(MT)]" \
	| sed 's/INS:NOVEL/INS/g' \
        | sed 's/DUP:INT/INS/g' \
        | sed 's/DUP:TANDEM/INS/g' \
	| grep 'SVTYPE=DEL;\|SVTYPE=INS;') \
        | bgzip -c > $OUT

	tabix $OUT
	printf "DONE\n"

	# Separate Insertions and Deletions, replace end position of insertions with starting position
	printf "Separate Insertions and deletions...."
	cat <(zcat $OUT | grep "^#") <(zcat $OUT | grep -vE "^#" | grep "SVTYPE=DEL;") | bgzip -c > $DEL
	cat <(zcat $OUT | grep "^#") <(zcat $OUT | grep -vE "^#" | grep "SVTYPE=INS;")  | bgzip -c > $INS

	# Index insertions and deletions
	tabix $DEL
	tabix $INS
	printf "DONE!\n"

	# ================================= Run TRUVARI for different qc values ============================================
	printf "Filter variants by quality score and call TRUVARI....\n"

	TMP="${FOLDER}/tmp"
	SUMMARY_INS="${FOLDER}/summary-insertions.tsv"
	SUMMARY_DEL="${FOLDER}/summary-deletions.tsv"
	touch "${SUMMARY_INS}"
	touch "${SUMMARY_DEL}"

	for ((qc=$2;qc<=$3;qc+=$4))
	do
		printf "\tProcessing score ${qc} of $3\n"

		# process deletions
		TMPQC="${TMP}/del_qc${qc}"
		mkdir -p "${TMPQC}"
		output="${TMPQC}/${qc}.vcf.gz"

		# filter by qc score and store file in tmp folder
		cat <(zcat $DEL | grep "^#") \
			<(zcat $DEL | grep -vE "^#" \
			| awk -v i="$qc" '$6>=i{print $0}') \
			| bgzip -c > "${output}"
		tabix "${output}"

		truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/deletions.vcf.gz \
			-c "$output" --includebed \
			/mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
			--passonly \
			-s 50 \
			-r 1000 \
			-p 0.00 \
			-o "${TMPQC}/truvari" &>/dev/null

		# extract precision and recall
		REC="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "recall" | cut -d':' -f2 | sed 's/,//g')"
		PREC="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "precision" | cut -d':' -f2 | sed 's/,//g')"
		F1="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "f1" | cut -d ':' -f 2 | sed 's/,//g')"
		printf "${qc}\t${PREC}\t${REC}\t${F1}\n" >> "${SUMMARY_DEL}"

		# process insertions
		TMPQC="${TMP}/ins_qc${qc}"
		mkdir -p "${TMPQC}"
		output="${TMPQC}/${qc}.vcf.gz"

		cat <(zcat $INS | grep "^#") \
			<(zcat $INS | grep -vE "^#" \
			| awk -v i="$qc" '$6>=i{print $0}') \
			| bgzip -c > "${output}"

		tabix "${output}"

		# call truvari on file

		truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/insertions.vcf.gz \
			-c "$output" --includebed \
			/mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
			--passonly \
			-s 50 \
			-r 1000 \
			-p 0.00 \
			-o "${TMPQC}/truvari" &>/dev/null

		# extract precision and recall
		REC="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "recall" | cut -d':' -f2 | sed 's/,//g')"
		PREC="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "precision" | cut -d':' -f2 | sed 's/,//g')"
		F1="$(cat "${TMPQC}/truvari/summary.txt" | grep -w "f1" | cut -d ':' -f 2 | sed 's/,//g')"
		printf "${qc}\t${PREC}\t${REC}\t${F1}\n" >> "${SUMMARY_INS}"
	done

	printf "DONE!\n"

	printf "Create plot...."
	$SCRIPTPATH/PrecRecplot.py --table $SUMMARY_INS --title Insertions
	$SCRIPTPATH/PrecRecplot.py --table $SUMMARY_DEL --title Deletions
	mv Deletions.png $FOLDER
	mv Insertions.png $FOLDER
	printf "DONE!\n"

	#printf "Cleanup...."
	# todo
	## delete tmp folder
	## maybe also delete extracted deletions and insertions
	#printf "DONE!\n"
fi
