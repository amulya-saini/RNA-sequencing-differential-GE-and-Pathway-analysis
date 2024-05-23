# Primary Sclerosing Cholangitis: From RNA Sequencing to Pathway Analysis

This project is a deep dive into the genomic landscapes governing arterial impact in Primary Sclerosing Cholangitis (PSC). 
Using RNA sequencing data derived from hepatic artery cells in liver transplant recipients. 
SRA Accession numbers:
"ERR9922199", "ERR9922200", "ERR9922201", "ERR9922202", "ERR9922203", "ERR9922204", "ERR9922205", "ERR9922206"

The study employs a robust analysis workflow encompassing tools like SRA-toolkit,fastqc, trimmomatic, Hisat2, subread, DESeq2, g:Profiler and KEGG pathway analysis. 
Results illuminate 297 differentially expressed genes, unraveling the molecular intricacies of PSC. 
The methodology involves rigorous preprocessing, quality assessment, and alignment steps, ensuring the reliability of downstream analyses. 
#### Noteworthy findings include 52 upregulated and 245 downregulated genes, providing critical insights into the genetic factors at play in PSC.


The code.txt file contains generalized commands for analyzing the data on a Unix system. However, I have also performed the analysis using Snakemake. Below are the instructions on how to do it. However, remember that the code should be tweaked to suit your analysis based on your data and specific requirements.

## What is Snakemake?

Snakemake is like our go-to tool for automating bioinformatics workflows. It's kind of like having a personal assistant for data analysis. Since we have to run a bunch of different programs in a certain order to process our data. Snakemake takes care of all that for us. We just need to tell it what tools we need, what data they work on, and how they fit together, and it handles the rest. It's super handy when we have a ton of files to process or when we're working on a big project that needs to run on a powerful computer. Basically, Snakemake helps us write down our analysis in a simple way, so we can focus on the science instead of worrying about the technical stuff.


## Snakemake Workflow for RNA-Seq Analysis

Below is a Snakemake workflow for RNA-Seq analysis. The workflow automates the processing of raw sequence data, including quality control, trimming, alignment, and quantification, using tools such as FastQC, Trimmomatic, Hisat2, and FeatureCounts. It is designed to run on a High-Performance Computing (HPC) system.

## Prerequisites

- Access to an HPC system
- Basic knowledge of HPC job submission (e.g., SLURM)
- Modules for required software tools should be available on the HPC system or installed manually (if applicable).

## Directory Structure

```plaintext
rna_seq_project/
├── data/
│   └── SRA/  # For raw SRA files
├── results/
│   ├── fastqc/
│   ├── trimmed/
│   ├── alignment/
│   ├── bam/
│   ├── sorted_bam/
│   ├── quantification/
├── logs/
├── scripts/
├── envs/
│   ├── sra.yaml
│   ├── fastqc.yaml
│   ├── multiqc.yaml
│   ├── trimmomatic.yaml
│   ├── hisat2.yaml
│   ├── samtools.yaml
│   ├── featurecounts.yaml
└── Snakefile
```
## Explanation

- `data/`: Directory for storing input data, including raw SRA files.
- `results/`: Directory for storing output files generated during the analysis.
- `logs/`: Directory for storing log files generated during the workflow execution.
- `scripts/`: Directory for storing helper scripts (if any).
- `envs/`: Directory containing Conda environment YAML files for different tools.
- `Snakefile`: Snakemake workflow definition file.

## Workflow Steps

1. **Create Conda Environments**: Set up Conda environments for each tool using the provided YAML files in the `envs/` directory.

2. **Initialize Conda Environments**:

```bash
conda env create -f envs/sra.yaml
conda env create -f envs/fastqc.yaml
conda env create -f envs/multiqc.yaml
conda env create -f envs/trimmomatic.yaml
conda env create -f envs/hisat2.yaml
conda env create -f envs/samtools.yaml
conda env create -f envs/featurecounts.yaml
```

3. **Run Snakemake Workflow**:

```bash
snakemake --use-conda --cores <num_cores>
```

Replace <num_cores> with the number of CPU cores you want Snakemake to utilize.

## Editing the Workflow
- Module Installation: Ensure that the required software tools are available on the HPC system or install them manually using Conda.
- Adjustments: Edit the Snakefile and YAML files according to your specific data and biological question.

## File Formats
- Conda Environment YAML Files: YAML files for each tool (e.g., fastqc.yaml, trimmomatic.yaml) should be saved in the envs/ directory.
- Snakefile: The Snakemake workflow definition file should be saved as Snakefile in the project root directory.
