

# 🧬 Tumor-Variant Calling Pipeline 
This pipeline performs tumor-variant calling using **GATK Mutect2** and **annotation using VEP (Variant Effect 

## 📂 Project Structure

- `Data/` → Raw sequencing data (FASTQ files)/NCBI/SRA
- `Scripts/` → Scripts for variant calling and annotation/Bash
- `Reports/` → Quality control (FastQC) and variant annotation (VEP) reports
- `Results/` → Processed files like BAM/SAM, VCF
- `Images/` → Visualizations of results


## 🛠️ Running the Pipeline

Below are the steps and commands to run the pipeline.

### 1. **Read Quality Control with FastQC**

To assess the qualityr FASTQ files, use **FastQC**:

```bash
fastqc data/SRR23020541_1_fastqc.html.fastq -o reports/fastqc/
``
2. Reference Genome Indexing with BWA
Before aligning the reads, We  need to index the reference genome (Homo_sapiens.GRCh38.dna.primary_assembly.fa) using BWA:

bwa index reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#This step creates the index files for reference genome (e.g., .amb, .ann, .bwt, .pac, .sa files).

3.. Read Alignment with BWA-MEM

bwa mem -t 4 reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
data/SRR23020541_1.fastq data/SRR23020541_2.fastq > results/SRR23020541_aligned.sam
After this exucation,  the indexed reference genome and paired-end FASTQ files to generate a SAM file. 

4. Convert SAM to BAM and Sort with Samtools:

By using samtools converted sam file to Bam file and sorted

samtools view -bS results/SRR23020541_aligned.sam > results/SRR23020541_aligned.bam

samtools sort results/SRR23020541_aligned.bam -o results/SRR23020541_sorted.bam

samtools index results/SRR23020541_sorted.bam

Step 1: Convert SAM to BAM
Step 2: Sort the BAM file by coordinates to prepare it for variant calling.
Step 3: Index the sorted BAM file for fast access during variant calling.

5. Variant Calling with GATK Mutect2/without normal sample

gatk Mutect2 -R reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I results/SRR23020541_sorted.bam -O results/SRR23020541_variants.vcf
 
vep -i results/SRR23020541_variants.vcf --cache --output_file reports/SRR23020541_vep_annotation.txt


#####📊 Results & Analysis

### 1. **Quality Control (FastQC)**

The quality of the raw sequencing reads was evaluated using **FastQC**. The FastQC reports provide detailed metrics such as base quality, GC content, sequence duplication levels, and adapter contamination. The **FastQC** results for the raw and trimmed FASTQ files are available as HTML reports and PNG images.

#### Raw Read QC:
- **FastQC Report (Raw)**: [SRR23020541_1_fastqc.html](reports/fastqc/SRR23020541_1_fastqc.html)
- **FastQC PNG (Raw)**: ![Raw FastQC](images/fastqc/SRR23020541_1_fastqc.png)

#### Trimmed Read QC:
After trimming the raw reads to remove low-quality bases and adapter sequences, another round of **FastQC** was performed on the trimmed reads to ensure the quality of the processed data.

- **FastQC Report (Trimmed)**: [SRR23020541_2_fastqc.html](reports/fastqc/SRR23020541_2_fastqc.html)
- **FastQC PNG (Trimmed)**: ![Trimmed FastQC](images/fastqc/SRR23020541_2_fastqc.png)

