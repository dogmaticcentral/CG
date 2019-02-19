awk -F '\t' '{print $8}' VHL_per_cancer_breakdown.tsv  | grep -v ENS | grep -v fs | grep -v pathogenic| \
 grep -v change | grep -v X | grep -v '*' | grep -v del |  grep -v '?' | grep -v PENY | grep -v RDAG | grep -v AGRP | \
 awk '{printf "%s   %s\n",  substr($1,0,1),  substr($1,2, length($1)-2)}' | sort -gk2 | uniq > tmp
