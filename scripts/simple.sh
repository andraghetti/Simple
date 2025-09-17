#!/bin/bash
#pipeline for mapping EMS mutants
#bug report to Lorenzo Andraghetti andraghetti.l@gmail.com

#mapping w/ BWA
echo "=========================================="
echo "STEP 0/8: SETTING UP RUN"
echo "=========================================="

# Create timestamped run directory
RUN_DIR="runs/run-$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_DIR"/{output,refs,fastq,scripts}

# Copy input files to run directory
cp -r fastq/* "$RUN_DIR/fastq/" 2>/dev/null || true

# Copy pipeline scripts to run directory for reproducibility
echo "Copying pipeline scripts to run directory for reproducibility..."
cp /app/scripts/simple.sh "$RUN_DIR/scripts/"
cp /app/scripts/simple_variables.sh "$RUN_DIR/scripts/"
cp /app/scripts/data_base.txt "$RUN_DIR/scripts/"
cp /app/scripts/analysis3.R "$RUN_DIR/scripts/"
echo "Scripts copied to $RUN_DIR/scripts/"

# Change to run directory for all operations
cd "$RUN_DIR"

{ # the entire script stdout and error will be displayed and redirected to log.txt
#reading variables
source ./scripts/simple_variables.sh

# Set Java memory allocation from configuration
JAVA_MEMORY=$java_memory
echo "Using Java memory allocation: $JAVA_MEMORY"

# Set CPU cores - auto-detect if set to "auto"
if [ "$cpu_cores" = "auto" ]; then
    CPU_CORES=$(nproc)
    echo "Auto-detected CPU cores: $CPU_CORES"
else
    CPU_CORES=$cpu_cores
    echo "Using specified CPU cores: $CPU_CORES"
fi

# Check for BWA - should be installed system-wide
if ! command -v bwa &> /dev/null; then
    echo "Error: bwa not found. This should have been installed during Docker build."
    exit 1
fi

# Check for samtools - should be installed system-wide
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. This should have been installed during Docker build."
    exit 1
fi


echo "Downloading & creating fasta file for $my_species"
fasta_link=`awk -v var="$my_species" 'match($1, var) {print $2}' ./scripts/data_base.txt`

if ! [ -f ./refs/$my_species.fa ]; then
  curl -o ./refs/$my_species.fa.gz $fasta_link
  gzip -d ./refs/$my_species.fa.gz
fi

if [ $my_species = "caenorhabditis_elegans" ]; then
    echo "Replacing roman numerals in numbers for $my_species"
    awk '{gsub(/^[\>][I]$/, ">1", $1); gsub(/^[\>][I][I]$/, ">2", $1); gsub(/^[\>][I][I][I]$/, ">3", $1); gsub(/^[\>][I][V]$/, ">4", $1); gsub(/^[\>][V]$/, ">5", $1);print}' ./refs/$my_species.fa > ./refs/$my_species.1.fa
fi

if [ -f ./refs/$my_species.1.fa ]; then 
    mv ./refs/$my_species.1.fa ./refs/$my_species.fa
fi

#head -n -1 (neg value) doesn't work on mac
awk '/[Ss]caffold/ || /[Cc]ontig/ {exit} {print}' ./refs/$my_species.fa > ./refs/$my_species.chrs.fa
fa=./refs/$my_species.chrs.fa

#choosing new release; 31 didn't work
echo "Downloading & creating knownsnps file for $my_species"
knownsnps_link=`awk -v var="$my_species" 'match($1, var) {print $3}' ./scripts/data_base.txt`
if ! [ -f ./refs/$my_species.vcf ]; then
  curl -o ./refs/$my_species.vcf.gz -Lk $knownsnps_link
  gzip -d ./refs/$my_species.vcf.gz
fi

#snpEff "link"
snpEff_link=`awk -v var="$my_species" 'match($1, var) {print $4}' ./scripts/data_base.txt`


#reference input files that are necessary to run the prograns
knownsnps=./refs/$my_species.vcf
#ftp://ftp.ensemblgenomes.org/pub/plants/release-31/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz
snpEffDB=$snpEff_link #paste the snpEff annotated genome name

echo "Creating .fai file for $my_species"
#creating .fai file
samtools faidx $fa
echo "Creating bwa index files for $my_species"
bwa index -p $my_species.chrs.fa -a is $fa
mv $my_species.chrs.* refs/

echo "Generating dict file for GATK"
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar CreateSequenceDictionary -R $fa -O refs/$my_species.chrs.dict


#making sure that all ref files where loaded and/or created; this is a good control and if the output is "something went wrong", it is usually picard failing due to a problem with java
# Expected files: 1 .fai + 5 BWA index files + 1 .dict = 7 files total
# Some genomes have additional files (ALT contigs, etc.) so more files is normal
a=`ls -l refs/ | wc -l`
if [ $a -ge 6 ]; then
	echo "$(tput setaf 2)refs loaded properly ($a files) $(tput setaf 7)"
else
	echo "$(tput setaf 1)something went wrong - expected at least 6 files, found $a $(tput setaf 7)"
	echo "Files in refs/:"
	ls -la refs/
fi


#mapping w/ BWA
echo "=========================================="
echo "STEP 1/8: BWA MAPPING"
echo "=========================================="
echo "This step maps reads to the reference genome..."
echo "Running BWA mapping with $CPU_CORES cores..."
bwa mem -t $CPU_CORES -M $fa ${mut_files[*]} > output/$mut.sam &
bwa mem -t $CPU_CORES -M $fa ${wt_files[*]} > output/$wt.sam
wait
echo "BWA mapping completed!"


#due to old samtools version this step is probably necessary
echo "Converting SAM to BAM format..."
samtools view -@ $CPU_CORES -bSh output/$mut.sam > output/$mut.bam &
samtools view -@ $CPU_CORES -bSh output/$wt.sam > output/$wt.bam
wait

rm -r output/*.sam

#this step is probably needed only when you have paired-end; in any case it should come before coordinate sorting (next step) on name-sorted files
echo "Fixing mate information..."
samtools fixmate -@ $CPU_CORES output/$mut.bam output/$mut.fix.bam &
samtools fixmate -@ $CPU_CORES output/$wt.bam output/$wt.fix.bam
wait

#sort by coordinates
echo "=========================================="
echo "STEP 2/8: SORTING BAM FILES"
echo "=========================================="
echo "This step sorts the BAM files by genomic coordinates..."
echo "Sorting BAM files..."
samtools sort -@ $CPU_CORES -o output/$mut.sort.bam output/$mut.fix.bam &
samtools sort -@ $CPU_CORES -o output/$wt.sort.bam output/$wt.fix.bam
wait
echo "BAM sorting completed!"

# Check if input BAM files exist before running MarkDuplicates
if [ ! -f output/$mut.sort.bam ]; then
    echo "ERROR: output/$mut.sort.bam not found!"
    exit 1
fi
if [ ! -f output/$wt.sort.bam ]; then
    echo "ERROR: output/$wt.sort.bam not found!"
    exit 1
fi

echo "=========================================="
echo "STEP 3/8: MARKING DUPLICATES"
echo "=========================================="
echo "This step identifies and marks PCR duplicates..."
echo "WARNING: This step can take a long time for large datasets!"
echo "Running MarkDuplicates for both samples..."
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar MarkDuplicates I=output/$mut.sort.bam O=output/$mut.sort.md.bam METRICS_FILE=output/$mut.matrics.txt ASSUME_SORTED=true &
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar MarkDuplicates I=output/$wt.sort.bam O=output/$wt.sort.md.bam METRICS_FILE=output/$wt.matrics.txt ASSUME_SORTED=true
wait

# Check if MarkDuplicates completed successfully
if [ ! -f output/$mut.sort.md.bam ]; then
    echo "ERROR: MarkDuplicates failed for mutant sample!"
    exit 1
fi
if [ ! -f output/$wt.sort.md.bam ]; then
    echo "ERROR: MarkDuplicates failed for wild-type sample!"
    exit 1
fi
echo "MarkDuplicates completed successfully for both samples."

#this part is just to add header for further gatk tools
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=output/$mut.sort.md.bam O=output/$mut.sort.md.rg.bam RGLB=$mut RGPL=illumina RGSM=$mut RGPU=run1 SORT_ORDER=coordinate &
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=output/$wt.sort.md.bam O=output/$wt.sort.md.rg.bam RGLB=$wt RGPL=illumina RGSM=$wt RGPU=run1 SORT_ORDER=coordinate
wait

$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar BuildBamIndex INPUT=output/$mut.sort.md.rg.bam OUTPUT=output/$mut.sort.md.rg.bai &
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/picard.jar BuildBamIndex INPUT=output/$wt.sort.md.rg.bam OUTPUT=output/$wt.sort.md.rg.bai
wait


#Variant calling using GATK HC extra parameters
echo "=========================================="
echo "STEP 5/8: VARIANT CALLING (GATK HaplotypeCaller)"
echo "=========================================="
echo "This step identifies genetic variants..."
echo "WARNING: This step can take a very long time (30+ minutes)!"
echo "Running GATK HaplotypeCaller with $CPU_CORES cores..."
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $fa -I output/$mut.sort.md.rg.bam -I output/$wt.sort.md.rg.bam -o output/$line.hc.vcf -minReadsPerAlignStart 7 -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES -writeFullFormat -stand_call_conf 10 -nct $CPU_CORES -variant_index_type LINEAR -variant_index_parameter 128000 -allowPotentiallyMisencodedQuals

# Check if HaplotypeCaller completed successfully
if [ ! -f output/$line.hc.vcf ]; then
    echo "ERROR: GATK HaplotypeCaller failed!"
    exit 1
fi
echo "GATK HaplotypeCaller completed successfully."

############preparing for R#########################
#Exclude indels from a VCF
#$java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar -R $fa -T SelectVariants --variant output/$line.hc.vcf -o output/$line.selvars.vcf --selectTypeToInclude SNP

echo "Converting VCF to table..."
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/GenomeAnalysisTK.jar -R $fa -T VariantsToTable -V output/$line.hc.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -o output/$line.table


####################################################################################################################################################
########################################now let's find the best candidates##########################################################################

#snpEff - Fix genome configuration
echo "=========================================="
echo "STEP 6/8: VARIANT ANNOTATION (snpEff)"
echo "=========================================="
echo "This step annotates variants with functional information..."
echo "WARNING: This step can take a long time (10+ minutes)!"
echo "Running snpEff annotation..."
$java -Xmx$JAVA_MEMORY -jar /usr/local/bin/snpEff.jar -c /usr/local/snpEff/snpEff.config -v $my_species output/$line.hc.vcf > output/$line.se.vcf 2> output/snpEff_error.log

# Check if snpEff succeeded
if [ ! -f "output/$line.se.vcf" ] || [ ! -s "output/$line.se.vcf" ]; then
    echo "WARNING: snpEff failed. Checking error log..."
    cat output/snpEff_error.log
    echo "Using original VCF for downstream analysis..."
    cp output/$line.hc.vcf output/$line.se.vcf
else
    echo "snpEff completed successfully!"
fi

###%%%%%%% JEN %%%%%%%%%%
#and finally, get only the SNPs that are ref/ref or ref/alt in the wt bulk and alt/alt in the mut bulk for recessive mutations
#for the case of dominant mutations should be ref/ref in the wt bulk and ref/alt or alt/alt in the mutant bulk
#column 10 is mutant bulk
#column 11 is WT bulk

echo "=========================================="
echo "STEP 7/8: FILTERING CANDIDATES"
echo "=========================================="
echo "This step filters variants based on mutation type and genotype patterns..."
echo "Filtering candidates based on mutation type: $mutation"

# Break down the long regex pattern to avoid "exceeds implementation size limit" error
# Define the effect patterns separately
EFFECT_PATTERNS="splice_acceptor_variant|splice_donor_variant|splice_region_variant|stop_lost|start_lost|stop_gained|missense_variant|coding_sequence_variant|inframe_insertion|disruptive_inframe_insertion|inframe_deletion|disruptive_inframe_deletion|exon_variant|exon_loss_variant|duplication|inversion|frameshift_variant|feature_ablation|gene_fusion|bidirectional_gene_fusion|rearranged_at_DNA_level|miRNA|initiator_codon_variant|start_retained"

if [ $mutation = "recessive" ]; then
	echo "Processing recessive mutations..."
	grep -v '^##' output/$line.se.vcf | awk -v pattern="$EFFECT_PATTERNS" 'BEGIN{FS=" "; OFS=" "} $1~/#CHROM/ || $10~/^1\/1/ && ($11~/^1\/0/ || $11~/^0\/0/ || $11~/^0\/1/) && $1~/^[0-9X]*$/ && $0~pattern {$3=$7=""; print $0}' | sed 's/  */ /g' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2] || $1~/#CHROM/) print $0}' > output/$line.cands2.txt
else
	echo "Processing dominant mutations..."
	grep -v '^##' output/$line.se.vcf | awk -v pattern="$EFFECT_PATTERNS" 'BEGIN{FS=" "; OFS=" "} $1~/#CHROM/ || ($10~/^0\/1/ || $10~/^1\/0/ || $10~/^1\/1/) && $11~/^0\/0/ && $1~/^[0-9X]*$/ && $0~pattern {$3=$7=""; print $0}' | sed 's/  */ /g' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2] || $1~/#CHROM/) print $0}' > output/$line.cands2.txt
fi


#getting things a bit more organized and only the relevant data from cands3
echo "Processing candidate data formatting..."
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" "$mut.ref" "$mut.alt" "$wt.ref" "$wt.alt" > output/$line.cands44.txt
echo "Formatting candidate data (this may take a moment for large datasets)..."
awk 'BEGIN{OFS="\t"} NR>1 {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' output/$line.cands2.txt | awk '$0!~/\./ && (($10+$11)>4) && (($12+$13)>4)' >> output/$line.cands44.txt

# JEN changed file name below
sort -k1,1 -k2 -n output/$line.cands44.txt > output/$line.candidates.txt

#this command will make it ready to run w/ R to produce Manhatten plot
echo "Preparing data for Manhattan plot generation..."
echo "Processing variant table data (this may take a moment)..."
printf "%s\t" "CHR" "POS" "REF" "ALT" "mut_GT" "mut.ref" "mut.alt" "mut.DP" "mut.GQ" "wt.GT" "wt.ref" "wt.alt" "wt.DP" "wt.GQ" > output/$line.plot.txt; printf "\n" >> output/$line.plot.txt
echo "Filtering variants for plotting..."
awk '$1~/^[0-9X]*$/ && $5~/^[AGCT]/ && $9~/^[AGCT]/ && $0 !~ /NA/ && $2 !~ /\./ && $3 !~ /\./ {gsub(/\,/, "\t"); print}' output/$line.table | awk '$6+$11>0 && $8>3 && $13>3' >> output/$line.plot.txt

#and finally, just get rid of known snps
echo "Removing known SNPs from plot data..."
awk 'FNR==NR{a[$1$2$4$5];next};!($1$2$3$4 in a)' $knownsnps output/$line.plot.txt > output/$line.plot.no_known_snps.txt

#get the snps in SnpEff format
echo "Processing SnpEff annotation data..."
echo "Matching variants with SnpEff annotations (this may take a moment)..."
awk 'FNR==NR{a[$1$2];next};($1$2 in a)' output/$line.plot.no_known_snps.txt output/$line.se.vcf > output/$line.plot2.txt
echo "Formatting SnpEff data..."
awk '{$3=$7=""; print $0}' output/$line.plot2.txt | sed 's/  */ /g' > output/$line.plot3.txt
echo "Filtering valid variants..."
awk '$3!~/\./ && $4!~/\./' output/$line.plot3.txt > output/$line.plot33.txt
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" "mut.ref" "mut.alt" "wt.ref" "wt.alt" > output/$line.plot44.txt
echo "Final formatting of plot data..."
awk 'BEGIN{OFS="\t"} {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' output/$line.plot33.txt >> output/$line.plot44.txt

##JEN changed filename below
echo "Sorting final plot data..."
sort -k1,1 -k2 -n output/$line.plot44.txt > output/$line.allSNPs.txt

#print cands that originate from a non-ref nucleotide
echo "Processing alternative allele candidates..."
echo "WARNING: This step processes large VCF files and may take several minutes!"
echo "Filtering alternative allele variants..."

grep -v '^##' output/$line.se.vcf | awk -v pattern="$EFFECT_PATTERNS" 'BEGIN{FS=" "; OFS=" "} $1~/#CHROM/ || ($10~/^2\/2/ && $11!~/^2\/2/) && $1~/^[0-9X]*$/ && $0~pattern {$3=$7=""; print $0}' | sed 's/  */ /g' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2] || $1~/#CHROM/) print $0}' > output/$line.cands_alt2.txt
echo "Removing known SNPs from alternative candidates..."
awk 'FNR==NR{a[$1$2$4$5];next};!($1$2$3$4 in a) || $1~/#CHROM/' $knownsnps output/$line.cands_alt2.txt > output/$line.cands_alt3.txt

#getting things a bit more organized and only the relevant data from cands3
echo "Formatting alternative allele candidate data..."
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" "$mut.ref" "$mut.alt" "$wt.ref" "$wt.alt" > output/$line.cands_alt4.txt
echo "Processing alternative allele formatting (this may take a moment)..."
awk 'BEGIN{OFS="\t"} NR>1 {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' output/$line.cands_alt3.txt | awk '$0!~/\./ && (($10+$11)>4) && (($12+$13)>4)' >> output/$line.cands_alt4.txt


#JEN added the line argument below
echo "=========================================="
echo "STEP 8/8: R ANALYSIS & PLOTTING"
echo "=========================================="
echo "This step generates Manhattan plots and final analysis..."
echo "Running R analysis..."
R --slave -e "library(ggplot2); library(reshape2); cat('R packages loaded successfully\n')" || {
    echo "ERROR: R packages not available!"
    exit 1
}
Rscript ./scripts/analysis3.R $line ./output

echo "=========================================="
echo "=========================================="
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "=========================================="
echo "All 8 steps completed successfully!"
echo "Results are available in the output directory."
echo "CPU cores used: $CPU_CORES"
echo "Java memory allocated: $JAVA_MEMORY"
echo "=========================================="
echo "Results saved in run directory: $(pwd)"
echo "Final results in: $(pwd)/output/"
echo "Performance: Used $CPU_CORES CPU cores for maximum speed"
echo "Memory: Allocated $JAVA_MEMORY RAM for Java tools"
echo "=========================================="

} 2>&1 | tee ./output/log.txt #with the { at the beginning of the text will redirect all console output to a file and still be visible in the terminal 




