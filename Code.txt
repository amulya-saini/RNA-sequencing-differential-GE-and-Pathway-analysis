
**Project Overview: Investigating Genomic Dynamics in Primary Sclerosing Cholangitis**

This project is a deep dive into the genomic landscapes governing arterial impact in Primary Sclerosing Cholangitis (PSC). Using RNA sequencing data derived from hepatic artery cells in liver transplant recipients, the study employs a robust analysis workflow encompassing tools like SRA-toolkit, Hisat2, DESeq2, and g:Profiler. Results illuminate 297 differentially expressed genes, unraveling the molecular intricacies of PSC.

The methodology involves rigorous preprocessing, quality assessment, and alignment steps, ensuring the reliability of downstream analyses. Noteworthy findings include 52 upregulated and 245 downregulated genes, providing critical insights into the genetic factors at play in PSC

Unix Commands used:

###### download SRA files ###### SRA TOOLKIT

- load the module
module load sra-toolkit

-Download the data
prefetch ascension_numbers
fastq-dump ascension_numbers

###### Quality control ###### FASTQC , MULTIQC

- loading fastqc
module load fastqc
fastqc path_to_the_fastq_files

- combining all the fastqc reports using multiqc
- Install the module
module load python
pip install multiqc

- add the module to path (in .bashrc) 
- navigate to the folder in which all the fastqc reports are present using cd
- once in the folder run 
multiqc .

###### Trimming ###### Trimmomatic

- loading module
module load trimmomatic

##trimming using a sliding window and no specific minimum length

# Input directory containing the original forward and reverse fastq files
input_dir="all_samples/"

# Output directory to store trimmed fastq files
output_dir="sliding_trimmed/"


# Iterate over each paired-end file in the input directory
for forward_file in "$input_dir"/*_1.fastq; do
    # Extract the file name without the extension
    filename=$(basename -- "$forward_file")

    # Define the corresponding reverse file
    reverse_file="$input_dir/${filename%_1.fastq}_2.fastq"

    # Trim using Trimmomatic for paired-end data without adapter removal
    trimmomatic PE -phred33 \
    "$forward_file" "$reverse_file" \
    "$output_dir/${filename%_1.fastq}_forward_paired.fastq" "$output_dir/${filename%_1.fastq}_forward_unpaired.fastq" "$output_dir/${filename%_1.fastq}_reverse_paired.fastq" "$output_dir/${filename%_1.fastq}_reverse_unpaired.fastq" \
    LEADING:3 TRAILING:0 SLIDINGWINDOW:4:20 > trim_log.txt 2>&1
done

- running the fastqc and multiqc 

###### Alignment ###### Hisat2

- loading module
module load hisat2

# Specify paths
index_folder="/N/slate/amsaini/reference_genome/hisat2_index"
input_folder="/N/slate/amsaini/Precision_project/sliding_trimmed"
output_folder="/N/slate/amsaini/Precision_project/hisat2_alignment"

# List of input fastq files (assuming paired-end data)
input_files=(
    "ERR9922199_forward_paired.fastq ERR9922199_reverse_paired.fastq"
    "ERR9922200_forward_paired.fastq ERR9922200_reverse_paired.fastq"
    "ERR9922201_forward_paired.fastq ERR9922201_reverse_paired.fastq"
    "ERR9922202_forward_paired.fastq ERR9922202_reverse_paired.fastq"
    "ERR9922203_forward_paired.fastq ERR9922203_reverse_paired.fastq"
    "ERR9922204_forward_paired.fastq ERR9922204_reverse_paired.fastq"
    "ERR9922205_forward_paired.fastq ERR9922205_reverse_paired.fastq"
    "ERR9922206_forward_paired.fastq ERR9922206_reverse_paired.fastq"
)

# Loop through samples
for files_pair in "${input_files[@]}"; do
    # Split each pair into an array
    files=($files_pair)
    R1="${files[0]}"
    R2="${files[1]}"

    # Extract sample name without extension
    sample_name=$(basename "${R1}" "_forward_paired.fastq")

    # Specify output alignment file
    output_alignment="${output_folder}/${sample_name}.sam"

    # Run hisat2
    hisat2 -x "${index_folder}/index" -1 "${input_folder}/${R1}" -2 "${input_folder}/${R2}" -S "${output_alignment}"
done

###### Sam to Bam conversion ###### SAMTOOLS

- loading modules
module swap gcc/12.1.0 gcc/9.3.0
module load python
module load samtools

# Specify the directories
sam_files_dir="/N/slate/amsaini/Precision_project/alignment_sam"
bam_files_dir="/N/slate/amsaini/Precision_project/bam_files"
sorted_indexed_dir="/N/slate/amsaini/Precision_project/bam_files/sorted_indexed"

# Create the directories if they don't exist
mkdir -p "${bam_files_dir}"
mkdir -p "${sorted_indexed_dir}"

# Iterate over the range of file numbers
for i in {199..206}; do
    # Define the input and output file paths
    input_sam="${sam_files_dir}/ERR9922${i}.sam"
    output_bam="${bam_files_dir}/ERR9922${i}.bam"
    sorted_bam="${sorted_indexed_dir}/ERR9922${i}_sorted.bam"
    indexed_bam="${sorted_indexed_dir}/ERR9922${i}_sorted.bam.bai"

    # Convert SAM to BAM
    samtools view -bS -o "${output_bam}" "${input_sam}"

    # Sort BAM
    samtools sort -o "${sorted_bam}" "${output_bam}"

    # Index sorted BAM
    samtools index "${sorted_bam}"
done

###### Quantification ###### FEATURECOUNTS

- loading the module
module load subread

# Specify the path to the folder containing BAM files
bam_folder="/N/slate/amsaini/Precision_project/hisat2_alignment/bam_sort_index"

# Specify the paths to gene and miRNA annotation files
gene_gtf_file="/N/slate/amsaini/reference_genome/illumina_reference/Homo_sapiens/Annotation/Genes/genes.gtf"

# Specify the output folder
output_folder="/N/slate/amsaini/Precision_project/quantification"

# Run featureCounts for gene annotations with paired-end reads
gene_output_file="$output_folder/all_samples_gene_count.txt"
featureCounts -a "$gene_gtf_file" -o "$gene_output_file" -p "${bam_folder}"/*.bam

##################Moving the files to the local directory for further downstream analysis in R (DEseq2) #######################

I used DEseq2 to perform differentially expressed genes in the attached R code. 
- Drop columns ("Chr", "Start", "End", "Strand", "Length") from the counts' files obtained from feature counts
- Change the first column to row names
- Combine both the datasets without losing row/column values
- filter the data frame by summing the counts across each and dropping rows(genes) with less than 1 count.
- Create DEseq dataset from the matrix using the metadata from the study
- Run DEseq Analysis
- Sort the adjusted p-values in ascending order
- Genes with less than 0.05 adjusted p_value and a log2foldchange of 1 (negative and positive) are selected as Differentially expressed genes


Data used:

- I have used Human reference genome version 38 (hg38), both fasta file and gtf file were downloaded from the following website(illumina):
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz

- The RNA seq data was obtained from the following study
htps://www.ncbi.nlm.nih.gov/bioproject/PRJEB54582


