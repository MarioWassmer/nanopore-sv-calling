 cat <(cat $1 | grep "^#") <(cat $1 | grep -vE "^#" | awk '$7 !~ /PASS/{print $0}') > nopass.vcf
