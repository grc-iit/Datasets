# Genomics Data Files

Standard genomics data formats for sequencing, alignment, and variant analysis used across next-generation sequencing (NGS) workflows.

**Documentation:** [biopython.org](http://biopython.org/) | **Format Specs:** [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/)

## Available Datasets

### synthetic_genome_reference.fasta
- **Size:** 179 KB
- **Type:** Genome reference (FASTA)
- **Description:** Synthetic reference genome with 5 chromosomes, 180 Kbp total

### gene_sequences.fasta
- **Size:** 1.1 KB
- **Type:** Gene sequences (FASTA)
- **Description:** Annotated gene sequences

### variant_calls.h5
- **Size:** 1.0 MB
- **Type:** Genetic variants (HDF5)
- **Description:** 5,000 SNPs/indels across 100 samples with genotypes

### rnaseq_expression.h5
- **Size:** 7.3 MB
- **Type:** RNA-seq expression data (HDF5)
- **Description:** 15,000 genes × 48 samples with differential expression data

### toy_alignment.sam
- **Size:** 786 bytes
- **Type:** Sequence alignment (SAM format)
- **Description:** Example aligned sequencing reads

### sample_variants.vcf
- **Size:** 2.1 KB
- **Type:** Variant calls (VCF format)
- **Description:** SNPs and indels in standard VCF format

### illumina_reads_sample.fastq
- **Size:** 178 bytes
- **Type:** Raw reads (FASTQ format)
- **Description:** Illumina sequencing reads with quality scores

## Where to Find More Data

- **1000 Genomes Project:** [internationalgenome.org](http://www.internationalgenome.org/) - Human genetic variation
- **NCBI SRA:** [ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra) - Sequence Read Archive (PB-scale)
- **TCGA:** [cancer.gov/tcga](https://www.cancer.gov/tcga) - Cancer genomics (TB-scale)
- **gnomAD:** [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/) - Human variants (>700TB)
- **ENCODE:** [encodeproject.org](https://www.encodeproject.org/) - Functional genomics
- **UK Biobank:** [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/) - 500k genomes
- **EMBL-EBI:** [ebi.ac.uk](https://www.ebi.ac.uk/) - European bioinformatics data
