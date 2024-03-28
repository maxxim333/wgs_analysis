# wgs_analysis
I analysed my father's whole genome sequencing data and made cool graphs

First, from the .VCF file, I extracted only those that have at least one alternate allele, using bcftools:

####Extract only variants that are homozygous for alternate allele, heatero AND have at least one read\
bcftools view -i 'FILTER="PASS" && (GT="1/1" || GT="0/1") && (AC>0) ' /Volumes/T7/nebula_genomics/f9949e99-72f4-488d-b81c-f03defe4fadb.vcf > filtered_only_alternate.vcf\



#I downloaded the data from clinvar one by one for each chromosome from https://www.ncbi.nlm.nih.gov/variation/view and stored it in /Users/maxxim333/Desktop/projects/haplogroupsY/clin_var_1by1 I only selected SNPS with a certain consequence that interested me. I named the files with a specific pattern like "ass37_chr18", "ass37_chr19" etc...
\
#I want to join them in one file but the problem is that there is no information of what chromosome it came from, so I need to add a column to each file to indicate the chromosome it belongs to. I did it with a bash code:\

#!/bin/bash\
\
# Directory containing the files\
directory="/Users/maxxim333/Desktop/projects/haplogroupsY/clin_var_1by1"\
\
# Iterate over files in the directory\
for filename in "$directory"/*.txt; do\
    if [ -f "$filename" ]; then\
        # Extract 'X' value from the filename\
        x_value=$(basename "$filename" | sed 's/.*_chr\\([^.]*\\)\\.txt/\\1/')\
\
        # Add 'X' value as the first column\
        awk -v x="$x_value" -F'\\t' -v OFS='\\t' '\{print x, $0\}' "$filename" > "$filename.tmp"\
\
        # Replace the original file with the modified one\
        mv "$filename.tmp" "$filename"\
    fi\


#Then I joined them via \
clin_var_1by1 % cat *.txt > ass37_all.vcf\

#Just to make sure all the columns are tab separated, I did:\
awk -v FS='\\t' -v OFS='\\t' '\{print\}' ass37_all.vcf > ass37_all.vcf

Then I ran the python code
