

# Download refFlat from here
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/

sed -i 's/chrM/MT/g' Hsapiens_GRCh38_110/refFlat.txt 
sed -i 's/chr//g' Hsapiens_GRCh38_110/refFlat.txt

# Download NGScheckmate
https://github.com/parklab/NGSCheckMate/blob/master/SNP/SNP_GRCh38_hg38_wChr.bed
sed -i 's/chrM/MT/g' Hsapiens_GRCh38_110/SNP_GRCh38_hg38_noChr.bed
sed -i 's/chr//g' Hsapiens_GRCh38_110/SNP_GRCh38_hg38_noChr.bed


# Create rRNA interval list
grep rRNA Homo_sapiens.GRCh38.110.chr.gtf | awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5}' > Homo_sapiens.GRCh38.110.chr.rRNA.bed

sed -i 's/chrM/MT/g' Hsapiens_GRCh38_110/refFlat.txt > Hsapiens_GRCh38_110/refFlat_noChr.txt 

# Create a genome bed file 
sudo apt install bedops
 gtf2bed < Homo_sapiens.GRCh38.110.chr.gtf > Homo_sapiens.GRCh38.110.chr.bed
 
 
# Chromosome sizes
pip install pyfaidx
faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa -i chromsizes > Homo_sapiens.GRCh38.dna.primary_assembly.fa.sizes