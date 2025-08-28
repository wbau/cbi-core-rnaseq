[README.md](https://github.com/user-attachments/files/22025852/README.md)

# cbi-core-rnaseq

An analysis pipeline for **total RNA-seq** data (Illumina paired-end reads).  
It performs read quality control, adapter/quality trimming, rRNA removal, alignment, gene quantification, and downstream QC aggregation.

---

## Pipeline summary

The pipeline implements the following steps:

1. **Input parsing**
   - Accepts a CSV samplesheet with sample names, FASTQ files, and strandedness.
   - Supports replicate handling: replicates are automatically merged into single FASTQ pairs before processing.

2. **Quality control (raw reads)**
   - **FastQC** to evaluate raw FASTQ files.

3. **Filtering and trimming**
   - **BBDuk** for adapter and quality trimming.

4. **rRNA removal**
   - **RiboDetector** to identify and remove rRNA contamination.

5. **Quality control (filtered reads)**
   - **FastQC** on filtered FASTQ files.

6. **Alignment and quantification**
   - **STAR** for genome alignment and gene-level counts (`--quantMode GeneCounts`).
   - Produces both coordinate-sorted BAM files and gene count tables.

7. **QC metrics**
   - **Picard**: RNA-seq metrics and insert size metrics.
   - **Samtools**: alignment statistics.
   - **RSeQC**: library complexity, read distribution, duplication, insert size, inner distance, splice junction metrics, etc.

8. **Reporting**
   - **MultiQC** aggregates results from all tools into a single HTML report.
   - **R (DESeq2, edgeR)** optional downstream differential expression and normalization (VST/CPM).

---

## Tools used

- FastQC
- BBDuk (BBMap suite)
- RiboDetector
- STAR
- Picard
- Samtools
- RSeQC
- MultiQC
- R + Bioconductor: DESeq2, edgeR

---

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/cbi-core-rnaseq.git
cd cbi-core-rnaseq
```

Install [Nextflow](https://www.nextflow.io/) (≥21.10.0) and ensure you have a working **Conda** or **Mamba** setup.  
Each process runs in its own Conda environment as defined in the pipeline.

---

## Input samplesheet

The pipeline requires a CSV file with the following columns:

- `sample` – Sample ID (replicates must have the same ID; they will be merged).
- `fastq_1` – Path to R1 FASTQ file (gzip-compressed).
- `fastq_2` – Path to R2 FASTQ file.
- `strandedness` – Library strandedness (e.g., `auto`, `forward`, `reverse`).

### Example

```csv
sample,fastq_1,fastq_2,strandedness
SAMPLE1_Tn,/data/SAMPLE1_S1_R1.fastq.gz,/data/SAMPLE1_S1_R2.fastq.gz,auto
SAMPLE1_Tn,/data/SAMPLE1_S2_R1.fastq.gz,/data/SAMPLE1_S2_R2.fastq.gz,auto
SAMPLE2_Tc,/data/SAMPLE2_S1_R1.fastq.gz,/data/SAMPLE2_S1_R2.fastq.gz,auto
SAMPLE3_Tc,/data/SAMPLE3_S1_R1.fastq.gz,/data/SAMPLE3_S1_R2.fastq.gz,auto
SAMPLE3_Tc,/data/SAMPLE3_S2_R1.fastq.gz,/data/SAMPLE3_S2_R2.fastq.gz,auto
SAMPLE4_blood,/data/SAMPLE4_blood_R1.fastq.gz,/data/SAMPLE4_blood_R2.fastq.gz,auto
```

Here:
- `SAMPLE1_Tn` has **two replicates**, which will be merged before analysis.
- `SAMPLE2_Tc` has a single pair.
- `SAMPLE3_Tc` has **two replicates** merged into one sample.
- `SAMPLE4_blood` is a separate blood sample.

---

## Running the pipeline

Basic run:

```bash
nextflow run main.nf \
   --samplesheet input/samplesheet.csv \
   --genome GENOME_NAME \
   --genome_fasta genome.fasta \
   --gtf genome.gtf \
   --bed genome.bed12 \
   --rRNA rrna_databases/ \
   --outdir results/rnaseq_project/
```

Resume from the last successful step:

```bash
nextflow run main.nf -resume --samplesheet input/samplesheet.csv
```

---

## Outputs

- `Reports/` – FastQC, MultiQC, Picard, RSeQC summaries.
- `Filtering/` – Filtered reads from BBDuk and RiboDetector.
- `Alignment/STAR/` – Aligned BAM files.
- `Counts/` – Gene count matrices (`*.ReadsPerGene.out.tab`).
- `R/` – Normalized counts, DESeq2 results, optional plots.

The final **MultiQC report** combines QC metrics from all steps into one interactive HTML file.

---

## Citation

If you use this pipeline, please cite:

- STAR: Dobin et al., Bioinformatics (2013).
- DESeq2: Love et al., Genome Biology (2014).
- edgeR: Robinson et al., Bioinformatics (2010).
- Plus the tools listed in the *Tools used* section above.
