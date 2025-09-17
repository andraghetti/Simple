# Simple - A Pipeline for Mapping Point Mutations

[![Ubuntu Base](https://img.shields.io/badge/Ubuntu-20.04-orange?&logo=ubuntu)](https://ubuntu.com/)
[![Java Version](https://img.shields.io/badge/OpenJDK-1.8.0_452-%23ED8B00?&logo=openjdk)](https://openjdk.org/projects/jdk8/)

[![BWA Version](https://img.shields.io/badge/BWA-0.7.19-green?&logo=github)](https://github.com/lh3/bwa/releases/tag/v0.7.19)
[![Samtools Version](https://img.shields.io/badge/Samtools-1.22.1-blue?&logo=github)](https://www.htslib.org/doc/1.22/)
[![GATK Version](https://img.shields.io/badge/GATK-v3.7--0-red?&logo=github)](https://github.com/broadgsa/gatk/releases/tag/3.7)
[![Picard Version](https://img.shields.io/badge/Picard-2.27.5-yellow?&logo=github)](https://github.com/broadinstitute/picard/releases/tag/2.27.5)
[![snpEff Version](https://img.shields.io/badge/snpEff-4.3p-purple?&logo=github)](https://sourceforge.net/projects/snpeff/files/snpEff_v4_3p_core.zip)

Simple is a bioinformatics pipeline for mapping EMS-induced point mutations using bulk segregant analysis. This Dockerized version provides an easy-to-use, reproducible environment for running the analysis.


## Complete Pipeline Walkthrough

### Step 1: Prepare Your Data

> [!IMPORTANT]
> **File Naming Convention**
> 
> Your FastQ files **MUST** follow this exact naming pattern in the `fastq/` directory:

### Required File Names:
```
fastq/
├── mut.R1.fastq    # Mutant bulk, read 1
├── mut.R2.fastq    # Mutant bulk, read 2 (paired-end)
├── wt.R1.fastq     # Wild-type bulk, read 1
└── wt.R2.fastq     # Wild-type bulk, read 2 (paired-end)
```

### Alternative Naming (with line prefix):
```
fastq/
├── YOUR_LINE_NAME.mut.R1.fastq
├── YOUR_LINE_NAME.mut.R2.fastq
├── YOUR_LINE_NAME.wt.R1.fastq
└── YOUR_LINE_NAME.wt.R2.fastq
```

### ⚠️ Common Mistakes to Avoid:
- ❌ `mutant.fastq` - Wrong! Must be `mut.R1.fastq`
- ❌ `wildtype.fastq` - Wrong! Must be `wt.R1.fastq`
- ❌ `mut_1.fastq` - Wrong! Must be `mut.R1.fastq`
- ❌ `wt_1.fastq` - Wrong! Must be `wt.R1.fastq`
- ❌ Missing `.R1` or `.R2` - Wrong! Required for paired-end

### ✅ Correct Examples:
- `mut.R1.fastq`, `mut.R2.fastq`, `wt.R1.fastq`, `wt.R2.fastq`
- `root_mutant.mut.R1.fastq`, `root_mutant.mut.R2.fastq`, `root_mutant.wt.R1.fastq`, `root_mutant.wt.R2.fastq`

### Step 2: Check Available Species
Before configuring, check what species are available in `scripts/data_base.txt`:

```bash
# View available species (first column)
cat scripts/data_base.txt | grep -v "^#" | awk '{print $1}'
```

**Available species include:**
- `Arabidopsis_thaliana`
- `Oryza_sativa_Japonica` 
- `Zea_mays`
- `Solanum_lycopersicum`
- `Drosophila_melanogaster`
- `caenorhabditis_elegans`
- `danio_rerio`
- `Saccharomyces_cerevisiae`

### Step 3: Configure the Analysis
Edit `scripts/simple_variables.sh` to set your parameters:

```bash
# Change the line name (used as prefix for output files, default is EMS)
line=YOUR_LINE_NAME

# Set mutation type
mutation=recessive  # or "dominant"

# Set species (must match entries in scripts/data_base.txt in organism column)
my_species=Arabidopsis_thaliana

# CPU cores configuration
cpu_cores="auto"  # Use "auto" for all available cores, or specify a number (e.g., 4)

# Java memory allocation
java_memory="8g"  # Set RAM for Java tools (e.g., "8g", "16g", "4g")
```

### Step 4: Run the Pipeline

#### Option A: Use Pre-built Docker Image (Recommended)
```bash
# Pull the pre-built image from GitHub Container Registry
# No login required - images are publicly accessible
# Docker automatically detects your platform (AMD64/ARM64) and pulls the correct image
docker pull ghcr.io/andraghetti/Simple:latest

# Run the pipeline
docker run --rm \
  -v $(pwd)/fastq:/app/fastq:ro \
  -v $(pwd)/runs:/app/runs \
  ghcr.io/andraghetti/Simple:latest
```

#### Option B: Build Docker Image Locally
```bash
# Build with explicit platform specification (works with regular docker build)
docker build -t simple-pipeline .

# Run the pipeline
docker run --rm \
  -v $(pwd)/fastq:/app/fastq:ro \
  -v $(pwd)/runs:/app/runs \
  simple-pipeline
```

### Step 5: Check Results
After completion, check the `runs/` directory for your results:
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.candidates.txt` - Candidate mutations
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.allSNPs.txt` - All SNPs for plotting
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.Rplot_*.pdf` - Manhattan plots
- `runs/run-YYYYMMDD_HHMMSS/output/log.txt` - Complete execution log
- `runs/run-YYYYMMDD_HHMMSS/archive/` - All intermediate files
- `runs/run-YYYYMMDD_HHMMSS/refs/` - Reference files used
- `runs/run-YYYYMMDD_HHMMSS/fastq/` - Input FASTQ files

## Features

- **Multi-Architecture Support**: Optimized for both AMD64 and ARM64 (Apple Silicon)
- **Dockerized Environment**: All dependencies pre-installed and compiled
- **Self-Contained Runs**: Each run creates a timestamped directory with all data
- **Pre-compiled Tools**: Latest versions of BWA, samtools, Picard, and snpEff
- **Multi-Core Performance**: Automatically detects and uses all available CPU cores
- **High Memory Allocation**: 8GB RAM (configurable) allocated for Java tools
- **Public Docker Images**: No build time required, freely accessible from GitHub Container Registry

### Pipeline Steps
1. **Reference Preparation**: Downloads and indexes reference genome
2. **Read Alignment**: Maps reads using BWA with optimized parameters
3. **BAM Processing**: Sorts, fixes mates, and marks duplicates
4. **Variant Calling**: Uses GATK HaplotypeCaller for variant discovery
5. **Variant Annotation**: Annotates variants with snpEff
6. **Candidate Selection**: Identifies candidate genes based on mutation type

## Example Analysis

### Root Development Mutant Analysis
1. **Prepare files**: `root_mutant.mut.R1.fastq`, `root_mutant.mut.R2.fastq`, `root_mutant.wt.R1.fastq`, `root_mutant.wt.R2.fastq`
2. **Configure**: Set `line=root_mutant`, `mutation=recessive`, `my_species=Arabidopsis_thaliana`
3. **Run**: Execute the Docker pipeline
4. **Results**: Check `root_mutant.candidates.txt` for candidate genes affecting root development

## Supported Species

The pipeline supports various species listed in `scripts/data_base.txt`. Common examples:
- Arabidopsis_thaliana
- Oryza_sativa_Japonica
- Zea_mays
- And many more...

## Troubleshooting

### Common Issues

1. **Pipeline fails immediately**: Check file naming convention in `fastq/` directory
2. **Out of memory errors**: Reduce `java_memory` in `simple_variables.sh`
3. **Slow performance**: Ensure `cpu_cores="auto"` for maximum performance
4. **Missing reference**: Check `my_species` matches entries in `data_base.txt`

### Getting Help

- Check the log files in `runs/run-TIMESTAMP/output/`
- Verify file naming convention matches exactly
- Ensure Docker has sufficient resources allocated

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite the original Simple paper and acknowledge this Dockerized version.