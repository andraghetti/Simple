#!/bin/bash
#pipeline for mapping EMS mutants
#bug report to Lorenzo Andraghetti andraghetti.l@gmail.com

# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -l, --line LINE_NAME        Line name (default: EMS)"
    echo "  -x, --mutation TYPE         Mutation type: recessive or dominant (default: recessive)"
    echo "  -s, --species SPECIES       Species name (default: Arabidopsis_thaliana)"
    echo "  -c, --cpu-cores CORES       CPU cores: auto or number (default: auto)"
    echo "  -m, --memory MEMORY         Memory allocation (default: auto)"
    echo "  -h, --help                  Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --line root_mutant --mutation recessive --species Arabidopsis_thaliana --cpu-cores auto --memory 16g"
    exit 1
}

# Default values
line="EMS"
mutation="recessive"
my_species="Arabidopsis_thaliana"
cores="auto"
memory="auto"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -l|--line)
            line="$2"
            shift 2
            ;;
        -x|--mutation)
            mutation="$2"
            shift 2
            ;;
        -s|--species)
            my_species="$2"
            shift 2
            ;;
        -c|--cpu-cores)
            cores="$2"
            shift 2
            ;;
        -m|--memory)
            memory="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate mutation type
if [[ "$mutation" != "recessive" && "$mutation" != "dominant" ]]; then
    echo "Error: mutation type must be 'recessive' or 'dominant'"
    exit 1
fi

# Handle "auto" memory setting
if [ "$memory" = "auto" ]; then
    # Auto-detect maximum safe memory (leave 2GB for system)
    total_mem_gb=$(free -g | awk '/^Mem:/{print $2}')
    if [ "$total_mem_gb" -gt 4 ]; then
        # Leave 2GB for system, use the rest
        safe_mem=$((total_mem_gb - 2))
        memory="${safe_mem}g"
    else
        # For systems with 4GB or less, use 2GB
        memory="2g"
    fi
fi

# Handle "auto" CPU cores setting
if [ "$cores" = "auto" ]; then
    # Auto-detect number of CPU cores
    cores=$(nproc)
fi

# Set derived variables
mut="${line}_mut"
wt="${line}_wt"

# Find input files with flexible naming patterns
# Supports: mut.R1.fq.gz, mut.R1.fastq, mut.R1.fq, etc.
# Also supports: line_mut.R1.fq.gz, line_mut.R1.fastq, etc.

echo "=========================================="
echo "STEP 0/8: SETTING UP RUN"
echo "=========================================="
echo "Configuration:"
echo "  Line: $line"
echo "  Mutation type: $mutation"
echo "  Species: $my_species"
echo "  CPU cores: $cores"
echo "  Java memory: $memory"
echo "=========================================="

# Create timestamped run directory
RUN_DIR="runs/run-$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_DIR"/{output,refs,data,scripts}

echo "Copying pipeline scripts to run directory for reproducibility..."
cp scripts/simple.sh "$RUN_DIR/scripts/"
cp scripts/data_base.txt "$RUN_DIR/scripts/"
cp scripts/analysis3.R "$RUN_DIR/scripts/"
echo "Scripts copied to $RUN_DIR/scripts/"
echo "Copying input files to run directory..."
cp data/* "$RUN_DIR/data/"
echo "Input files copied to $RUN_DIR/data/"


# Change to run directory for all operations
cd "$RUN_DIR"

# Search for input files in the copied data directory
echo "Searching for input files..."
mut_files=()
wt_files=()

# Look for mutant files (must contain "mut" and "R1" or "R2")
for file in data/*mut*R1* data/*mut*R2*; do
    if [ -f "$file" ]; then
        mut_files+=("$file")
        echo "Found mutant file: $file"
    fi
done

# Look for wild-type files (must contain "wt" and "R1" or "R2")
for file in data/*wt*R1* data/*wt*R2*; do
    if [ -f "$file" ]; then
        wt_files+=("$file")
        echo "Found wild-type file: $file"
    fi
done

# Validate that we found files
if [ ${#mut_files[@]} -eq 0 ]; then
    echo "ERROR: No mutant files found!"
    echo "Expected files containing 'mut' and 'R1' or 'R2' in data/ directory"
    echo "Examples: data/mut.R1.fq.gz, data/mut.R1.fastq, data/line_mut.R1.fq.gz"
    exit 1
fi

if [ ${#wt_files[@]} -eq 0 ]; then
    echo "ERROR: No wild-type files found!"
    echo "Expected files containing 'wt' and 'R1' or 'R2' in data/ directory"
    echo "Examples: data/wt.R1.fq.gz, data/wt.R1.fastq, data/line_wt.R1.fq.gz"
    exit 1
fi

echo "Found ${#mut_files[@]} mutant files and ${#wt_files[@]} wild-type files"

{ # the entire script stdout and error will be displayed and redirected to log.txt

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
java -Xmx$memory -jar /usr/local/bin/picard.jar CreateSequenceDictionary -R $fa -O refs/$my_species.chrs.dict


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
echo "Running BWA mapping with $cores cores..."
bwa mem -t $cores -M $fa ${mut_files[*]} > output/$mut.sam &
bwa mem -t $cores -M $fa ${wt_files[*]} > output/$wt.sam
wait
echo "BWA mapping completed!"


#due to old samtools version this step is probably necessary
echo "Converting SAM to BAM format..."
samtools view -@ $cores -bSh output/$mut.sam > output/$mut.bam &
samtools view -@ $cores -bSh output/$wt.sam > output/$wt.bam
wait

rm -r output/*.sam

#this step is probably needed only when you have paired-end; in any case it should come before coordinate sorting (next step) on name-sorted files
echo "Fixing mate information..."
samtools fixmate -@ $cores output/$mut.bam output/$mut.fix.bam &
samtools fixmate -@ $cores output/$wt.bam output/$wt.fix.bam
wait

#sort by coordinates
echo "=========================================="
echo "STEP 2/8: SORTING BAM FILES"
echo "=========================================="
echo "This step sorts the BAM files by genomic coordinates..."
echo "Sorting BAM files..."
samtools sort -@ $cores -o output/$mut.sort.bam output/$mut.fix.bam &
samtools sort -@ $cores -o output/$wt.sort.bam output/$wt.fix.bam
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
java -Xmx$memory -jar /usr/local/bin/picard.jar MarkDuplicates I=output/$mut.sort.bam O=output/$mut.sort.md.bam METRICS_FILE=output/$mut.matrics.txt ASSUME_SORTED=true &
java -Xmx$memory -jar /usr/local/bin/picard.jar MarkDuplicates I=output/$wt.sort.bam O=output/$wt.sort.md.bam METRICS_FILE=output/$wt.matrics.txt ASSUME_SORTED=true
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
java -Xmx$memory -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=output/$mut.sort.md.bam O=output/$mut.sort.md.rg.bam RGLB=$mut RGPL=illumina RGSM=$mut RGPU=run1 SORT_ORDER=coordinate &
java -Xmx$memory -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=output/$wt.sort.md.bam O=output/$wt.sort.md.rg.bam RGLB=$wt RGPL=illumina RGSM=$wt RGPU=run1 SORT_ORDER=coordinate
wait

java -Xmx$memory -jar /usr/local/bin/picard.jar BuildBamIndex INPUT=output/$mut.sort.md.rg.bam OUTPUT=output/$mut.sort.md.rg.bai &
java -Xmx$memory -jar /usr/local/bin/picard.jar BuildBamIndex INPUT=output/$wt.sort.md.rg.bam OUTPUT=output/$wt.sort.md.rg.bai
wait


#Variant calling using GATK HC extra parameters
echo "=========================================="
echo "STEP 5/8: VARIANT CALLING (GATK HaplotypeCaller)"
echo "=========================================="
echo "This step identifies genetic variants..."
echo "WARNING: This step can take a very long time (30+ minutes)!"
echo "Running GATK HaplotypeCaller with $cores cores..."
java -Xmx$memory -jar /usr/local/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $fa -I output/$mut.sort.md.rg.bam -I output/$wt.sort.md.rg.bam -o output/$line.hc.vcf -minReadsPerAlignStart 7 -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES -writeFullFormat -stand_call_conf 10 -nct $cores -variant_index_type LINEAR -variant_index_parameter 128000 -allowPotentiallyMisencodedQuals

# Check if HaplotypeCaller completed successfully
if [ ! -f output/$line.hc.vcf ]; then
    echo "ERROR: GATK HaplotypeCaller failed!"
    exit 1
fi
echo "GATK HaplotypeCaller completed successfully."

############preparing for R#########################
#Exclude indels from a VCF
java -Xmx2g -jar /usr/local/bin/GenomeAnalysisTK.jar -R $fa -T SelectVariants --variant output/$line.hc.vcf -o output/$line.selvars.vcf --selectTypeToInclude SNP

echo "Converting VCF to table..."
java -Xmx$memory -jar /usr/local/bin/GenomeAnalysisTK.jar -R $fa -T VariantsToTable -V output/$line.hc.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -o output/$line.table


####################################################################################################################################################
########################################now let's find the best candidates##########################################################################

#snpEff - Fix genome configuration
echo "=========================================="
echo "STEP 6/8: VARIANT ANNOTATION (snpEff)"
echo "=========================================="
echo "This step annotates variants with functional information..."
echo "WARNING: This step can take a long time (10+ minutes)!"
echo "Running snpEff annotation..."
java -Xmx$memory -jar /usr/local/bin/snpEff.jar -c /usr/local/snpEff/snpEff.config -v $my_species output/$line.hc.vcf > output/$line.se.vcf 2> output/snpEff_error.log

# Check if snpEff succeeded
if [ ! -f "output/$line.se.vcf" ] || [ ! -s "output/$line.se.vcf" ]; then
    echo "WARNING: snpEff failed. Checking error log..."
    cat output/snpEff_error.log
    echo "Using original VCF for downstream analysis..."
    cp output/$line.hc.vcf output/$line.se.vcf
else
    echo "snpEff completed successfully!"
fi

# Move snpEff summary files to output directory
echo "Organizing snpEff output files..."
if [ -f "snpEff_genes.txt" ]; then
    mv snpEff_genes.txt output/
fi
if [ -f "snpEff_summary.html" ]; then
    mv snpEff_summary.html output/
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
echo "Results saved in run directory: $(pwd)"
echo "Final results in: $(pwd)/output/"
echo "=========================================="

} 2>&1 | tee ./output/log.txt # this will redirect all console output to a file and still be visible in the terminal