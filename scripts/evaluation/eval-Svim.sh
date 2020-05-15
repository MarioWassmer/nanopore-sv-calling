#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "Usage: $0 sniffles_output.vcf[.gz] threshold(int in [1,100])"
	exit 1
else
	# Save input file in variable
	IN=$1
	THRESHOLD=$2

	# get absolute path of folder where the script has been started
	FOLDER="$(dirname $(readlink -f IN))/evaluation"
	mkdir -p $FOLDER

	# define output for survivor
	SURV_PREPARED=$FOLDER/$(echo $(basename $IN) | cut -f 1 -d '.')_survivor.vcf

	# define survivor output
	mkdir -p $FOLDER/survivor
	SURV_STAT=$FOLDER/survivor/survivor-stats.txt
	SURV_SUM=$FOLDER/survivor/survivor-summary.txt
	
	# define truvari prepared file
        TRUV_PREPARED=$FOLDER/$(echo $(basename $IN) | cut -f 1 -d '.')_ins-del.vcf.gz

	# define deletions and insertions output
	DEL=$FOLDER/deletions.vcf.gz
	INS=$FOLDER/insertions.vcf.gz

	# truvari output folder
	TRUV_DEL=$FOLDER/truvari-deletions
	TRUV_INS=$FOLDER/truvari-insertions

	# debug output
	#printf "Debug mode.. Printing all defined variables:\n"
	#printf "Input: $IN\n"
	#printf "Folder: $FOLDER\n"
	#printf "TRUV_PREPARED: $TRUV_PREPARED\n"
	#printf "DEL: $DEL\n"
	#printf "Ins: $INS\n"
	#printf "TRUV_DEL: $TRUV_DEL\n"


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

	# Create and output SURVIVOR statistics on the input file

	printf "Compute statistics...."
        cat <(zcat $IN | grep "^#") \
        <(zcat $IN | grep -vE "^#" \
	| egrep "^[(0-9)(X)(Y)(MT)]" \
	| awk -v threshold=$THRESHOLD '$6>=threshold {print $0}')\
       	> $SURV_PREPARED

	SURVIVOR stats $SURV_PREPARED 49 -1 -1 $SURV_STAT > $SURV_SUM
	printf "DONE!\n\n"
	
	cat $SURV_SUM

	printf "\n\n"


	# ================ Prepare for Truvari =================
	printf "Prepare input file for truvari evaluation...."
	# fix the header of the input file
	# Extract Insertions, Deletions and Duplications; Normalize Duplications to Insertions


	cat <(cat $SURV_PREPARED | grep "^#") \
        <(cat $SURV_PREPARED | grep -vE "^#" | grep -v "GL000\|hs37d5\|NC_\|MT" \
	| sed 's/INS:NOVEL/INS/g' \
        | sed 's/DUP:INT/INS/g' \
        | sed 's/DUP:TANDEM/INS/g' \
        | grep 'SVTYPE=DEL;\|SVTYPE=INS;' \
        | awk -v threshold=$THRESHOLD '$6>=threshold {print $0}') \
        | bgzip -c > $TRUV_PREPARED

	tabix $TRUV_PREPARED

	printf "DONE\n"

	# Now the file is fully prepared for truvari. Save it in current folder

	printf "Separate insertions  from deletions...."

	# Separate Insertions and Deletions, replace end position of insertions with starting position
        cat <(zcat $TRUV_PREPARED | grep "^#") \
		<(zcat $TRUV_PREPARED | grep -vE "^#" \
		| grep "SVTYPE=DEL;") \
		| bgzip -c > $DEL
	 
	tabix $DEL
         
	cat <(zcat $TRUV_PREPARED | grep "^#") \
		<(zcat $TRUV_PREPARED | grep -vE "^#" \
		| grep "SVTYPE=INS;") \
		| bgzip -c > $INS

	tabix $INS

	printf "DONE!\n"

	# ================ Run Truvari =========================
	
	# evaluate deletions
  	truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/deletions.vcf.gz \
		-c $DEL \
		--includebed /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
		-r 1000 \
		-p 0.00 \
		--passonly \
		-s 50 \
		-o $TRUV_DEL

	# evaluate deletions
	truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/insertions.vcf.gz \
        	-c $INS \
		--includebed /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
		-r 1000 \
		-p 0.00 \
		--passonly \
		-s 50 \
		-o $TRUV_INS

	# cleanup
	rm $SURV_PREPARED
fi
