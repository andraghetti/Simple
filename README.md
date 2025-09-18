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
> **Create a Working Directory**
> 
> **You must create a folder containing your FastQ files before running the pipeline.** The pipeline will create all output in a timestamped subdirectory within this folder.
> If you have cloned this repository, just use the `data` and `runs` directories.
> Otherwise create them inside a clean directory.

### Required File Names:
The pipeline automatically detects files based on flexible patterns. Files must contain:
- **"mut"** for mutant bulk samples
- **"wt"** for wild-type bulk samples  
- **"R1"** or **"R2"** to indicate read direction
- **Any extension**: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`, etc.

<details>

<summary>Click here to know more about the file structure and file naming</summary>

### ✅ Setup Example:
```bash
# 1. Create your analysis folder
mkdir my_mutant_analysis
cd my_mutant_analysis

# 2. Create data folder for your FastQ files
mkdir -p data runs

# 3. Copy your FastQ files to the data folder
# Your folder should now look like this:
ls -la data/
# mut.R1.fq.gz    mut.R2.fq.gz    wt.R1.fq.gz    wt.R2.fq.gz
```

### ✅ Supported Examples:
```
my_analysis/
├── data/
│   ├── mut.R1.fq.gz          # Mutant bulk, read 1 (compressed)
│   ├── mut.R2.fq.gz          # Mutant bulk, read 2 (compressed)
│   ├── wt.R1.fastq           # Wild-type bulk, read 1 (uncompressed)
│   └── wt.R2.fastq           # Wild-type bulk, read 2 (uncompressed)
└── runs/                     # Created automatically by pipeline
```

### ✅ Alternative Naming (with line prefix):
```
my_analysis/
├── data/
│   ├── root_mutant.mut.R1.fq.gz
│   ├── root_mutant.mut.R2.fq.gz
│   ├── root_mutant.wt.R1.fq.gz
│   └── root_mutant.wt.R2.fq.gz
└── runs/                     # Created automatically by pipeline
```

### ✅ More Examples:
```
my_analysis/
├── data/
│   ├── sample_mut.R1.fastq.gz
│   ├── sample_mut.R2.fastq.gz
│   ├── sample_wt.R1.fastq.gz
│   └── sample_wt.R2.fastq.gz
└── runs/                     # Created automatically by pipeline
```

### 📁 Repository Structure (for developers):
```
simple-fork/
├── data/                 # Example data folder (rename to your analysis name)
│   ├── mut.R1.fq.gz
│   ├── mut.R2.fq.gz
│   ├── wt.R1.fq.gz
│   └── wt.R2.fq.gz
├── scripts/              # Pipeline scripts
├── programs/             # Bioinformatics tools
└── Dockerfile
```

### ⚠️ Common Mistakes to Avoid:
- ❌ `mutant.fastq` - Missing "R1" or "R2" identifier
- ❌ `wildtype.fastq` - Missing "R1" or "R2" identifier  
- ❌ `mut_1.fastq` - Should be `mut.R1.fastq`
- ❌ `wt_1.fastq` - Should be `wt.R1.fastq`
- ❌ Files without "mut" or "wt" in the name
- ❌ Files without "R1" or "R2" in the name

### ✅ Correct Examples:
- `mut.R1.fq.gz`, `mut.R2.fq.gz`, `wt.R1.fq.gz`, `wt.R2.fq.gz`
- `root_mutant.mut.R1.fastq`, `root_mutant.mut.R2.fastq`, `root_mutant.wt.R1.fastq`, `root_mutant.wt.R2.fastq`
- `sample_mut.R1.fq`, `sample_mut.R2.fq`, `sample_wt.R1.fq`, `sample_wt.R2.fq`

</details>

### Step 2: Check Available Species
Before configuring, check what species are available:

```bash
# View available species using Docker
docker run --rm ghcr.io/andraghetti/simple cat /app/scripts/data_base.txt | grep -v "^#" | awk '{print $1}'
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

This will be required in the next step, so copy the value. It will be used in the `--species` attribute.
Example: `--species Drosophila_melanogaster`.

### Step 3: Configure the Analysis
The pipeline accepts command line arguments that override the default values. You can customize the analysis by passing parameters to the Docker container:

**Available Parameters:**
- `--line LINE_NAME`: Line name (default: EMS)
- `--mutation TYPE`: Mutation type: recessive or dominant (default: recessive)
- `--species SPECIES`: Species name (default: Arabidopsis_thaliana)
- `--cpu-cores CORES`: CPU cores: auto or number (default: auto)
- `--memory MEMORY`: Java memory allocation (default: auto)

### Step 4: Run the Pipeline

#### Option A: Use Pre-built Docker Image (Recommended)
```bash
# Pull the pre-built image from GitHub Container Registry
docker pull ghcr.io/andraghetti/simple:latest

# Run the pipeline with default settings
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/runs:/app/runs \
  ghcr.io/andraghetti/simple:latest
```

Alternative example with all the configurable options:

```bash
# Run the pipeline with custom parameters
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/runs:/app/runs \
  ghcr.io/andraghetti/simple:latest \
  --line root_mutant \
  --mutation recessive \
  --species Arabidopsis_thaliana \
  --cpu-cores 8 \
  --java-memory 16g
```

#### Option B: Build Docker Image Locally
```bash
# Build with explicit platform specification (works with regular docker build)
docker build -t simple-pipeline .

# Run the pipeline with default settings
docker run --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/runs:/app/runs \
  simple-pipeline
```

### Step 5: Check Results
After completion, check the `runs/` directory for your results:
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.candidates.txt` - Candidate mutations
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.allSNPs.txt` - All SNPs for plotting
- `runs/run-YYYYMMDD_HHMMSS/output/YOUR_LINE_NAME.Rplot_*.pdf` - Manhattan plots
- `runs/run-YYYYMMDD_HHMMSS/output/snpEff_summary.html` - snpEff annotation summary
- `runs/run-YYYYMMDD_HHMMSS/output/snpEff_genes.txt` - snpEff gene annotations
- `runs/run-YYYYMMDD_HHMMSS/output/log.txt` - Complete execution log


### Pipeline Steps
1. **Reference Preparation**: Downloads and indexes reference genome
2. **Read Alignment**: Maps reads using BWA with optimized parameters
3. **BAM Processing**: Sorts, fixes mates, and marks duplicates
4. **Variant Calling**: Uses GATK HaplotypeCaller for variant discovery
5. **Variant Annotation**: Annotates variants with snpEff
6. **Candidate Selection**: Identifies candidate genes based on mutation type

## Example Analysis

### Root Development Mutant Analysis
1. **Create analysis folder and add files**:
   ```bash
   # Create your analysis folder
   mkdir root_analysis
   cd root_analysis
   
   # Add your FastQ files (copy or download them here)
   # Files: root_mutant.mut.R1.fq.gz, root_mutant.mut.R2.fq.gz
   #        root_mutant.wt.R1.fq.gz, root_mutant.wt.R2.fq.gz
   ```

2. **Run the pipeline**:
   ```bash
   docker run --rm \
     -v $(pwd)/data:/app/data \
     -v $(pwd)/runs:/app/runs \
     ghcr.io/andraghetti/simple:latest \
     --line root_mutant \
     --mutation recessive \
     --species Arabidopsis_thaliana
   ```

3. **Check results**: Look in `runs/run-YYYYMMDD_HHMMSS/output/root_mutant.candidates.txt` for candidate genes affecting root development

## Troubleshooting

### Common Issues

1. **Pipeline fails immediately**: Check file naming convention in your analysis directory
   - Files must contain "mut" or "wt" AND "R1" or "R2"
   - Examples: `mut.R1.fq.gz`, `wt.R1.fastq`, `line_mut.R1.fq.gz`
2. **Out of memory errors**: Reduce `--memory` parameter (e.g., use `--memory 8g`)
3. **Slow performance**: Ensure `--cpu-cores auto` for maximum performance
4. **Missing reference**: Check `--species` matches entries in `data_base.txt`
5. **Invalid arguments**: Use `--help` to see available options
6. **No files found**: Verify files are in the current directory with correct naming pattern
7. **Permission errors**: Ensure Docker has read/write access to your analysis directory

### Getting Help

- Check the log files in `runs/run-TIMESTAMP/output/`
- Verify file naming convention matches exactly
- Ensure Docker has sufficient resources allocated

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite the original Simple paper and acknowledge this optimized version.
