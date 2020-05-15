#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "Usage: $0 sniffles_output.vcf[.gz] threshold]"
	exit 1
else
	# Save input file in variable
	IN=$1
	THRESHOLD=$2

	# get absolute path of folder where the script has been started
	FOLDER="$(dirname $(readlink -f IN))/evaluation"
	mkdir -p $FOLDER
	
	# define truvari prepared file
        TRUV_PREPARED=$FOLDER/$(echo $(basename $IN) | cut -f 1 -d '.')_ins-del.vcf.gz
        SURV_PREPARED=$FOLDER/$(echo $(basename $IN) | cut -f 1 -d '.')_survivor.vcf

	# define deletions and insertions output
	DEL=$FOLDER/deletions.vcf.gz
	INS=$FOLDER/insertions.vcf.gz

	# truvari output folder
	TRUV_DEL=$FOLDER/truvari-deletions
	TRUV_INS=$FOLDER/truvari-insertions

	# define survivor output
	mkdir -p $FOLDER/survivor
	SURV_STAT=$FOLDER/survivor/survivor-stats.txt
	SURV_SUM=$FOLDER/survivor/survivor-summary.txt

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

	cat <(zcat $IN | grep "^#" | sed 's/SUPTYPE,Number=A/SUPTYPE,Number=./g') \
	<(zcat $IN | grep -vE "^#" \
	| egrep "^[(0-9)(X)(Y)(MT)]" \
	| awk -v threshold="$THRESHOLD" '{n=split($8,w,";"); for(i=0;++i<=n;) if (w[i] ~ '/RE=/') split(w[i],depth," "); split(depth[1],re,"="); if(re[2]>=threshold) print $0}') \
      	> $SURV_PREPARED

	SURVIVOR stats $SURV_PREPARED 49 -1 -1 $SURV_STAT > $SURV_SUM
	printf "DONE!\n\n"

	cat $SURV_SUM

	printf "\n\n"

	# ================ Prepare for Truvari =================
	# fix the header of the input file
	# Extract Insertions, Deletions and Duplications; Normalize Duplications to Insertions


	cat <(cat $SURV_PREPARED | grep "^#") \
	<(cat $SURV_PREPARED | grep -vE "^#" \
	| grep "DUP\|INS\|DEL" \
	| sed "s/DUP/INS/g" \
	| sort -k1,1 -k2,2n) \
	| bgzip -c > $TRUV_PREPARED

	tabix $TRUV_PREPARED
	# Now the file is fully prepared for truvari. Save it in current folder

	 # Separate Insertions and Deletions, replace end position of insertions with starting position
         cat <(zcat $TRUV_PREPARED | grep "^#") \
		 <(zcat $TRUV_PREPARED | grep -vE "^#" \
		 | grep "SVTYPE=DEL;") \
		 | bgzip -c > $DEL
	 
	 tabix $DEL
         
	 cat <(zcat $TRUV_PREPARED | grep "^#") \
		 <(zcat $TRUV_PREPARED | grep -vE "^#" \
		 | grep "SVTYPE=INS;" \
		 | awk -v OFS='\t' '{gsub(/END=[0-9]+;/,"END="$2";",$8); print $0}') \
		 | bgzip -c > $INS

	 tabix $INS

	# ================ Run Truvari =========================
	
	# evaluate deletions
  	truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/deletions.vcf.gz \
		-c $DEL \
		--includebed /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
		-r 250 \
		-p 0.00 \
		--passonly \
		-s 50 \
		-o $TRUV_DEL

	# evaluate deletions
	truvari -b /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/evaluation/insertions.vcf.gz \
        	-c $INS \
		--includebed /mnt/SRV018/users/ahwassm1/Masterthesis/data/human/HG002_NA24385_son/hg002-tier1-v0.6/original/HG002_SVs_Tier1_v0.6.bed \
		-r 250 \
		-p 0.00 \
		--passonly \
		-s 50 \
		-o $TRUV_INS

	# =============== Output a summary for further evaluation ==============
	# todo


	# ============================================= Cleanup =====================================================
	rm $SURV_PREPARED
fi
