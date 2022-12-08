for f in *.tsv; do
    awk '{print $1","$4}' OFS='\t' "$f" > "/Users/henry/Documents/GitHub/Final/Gene_Counts_Modified/test.tsv" 
done