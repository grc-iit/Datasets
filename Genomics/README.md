# Genomics Data Files

Standard genomics data formats for sequencing, alignment, and variant analysis. These formats are used across next-generation sequencing (NGS) workflows, from raw reads to variant calling and genomic analysis.

## Available Datasets

### toy_alignment.sam
**Size:** 786 bytes
**Format:** SAM (Sequence Alignment/Map)
**Source:** samtools examples
**Type:** Aligned sequencing reads

Example alignment file showing mapped sequencing reads. SAM is the human-readable text format for storing biological sequences aligned to a reference genome.

**Use cases:** Alignment viewing, format learning, pipeline testing

### sample_variants.vcf
**Size:** 2.1 KB
**Format:** VCF (Variant Call Format)
**Type:** Genetic variants

Example variant calls showing SNPs, insertions, and deletions. VCF is the standard format for storing gene sequence variations.

**Use cases:** Variant analysis, genotype studies, population genetics

### illumina_reads_sample.fastq
**Size:** 178 bytes
**Format:** FASTQ
**Type:** Raw sequencing reads with quality scores

Raw sequencing reads from Illumina platform with Phred quality scores. FASTQ is the de facto standard for storing unaligned sequencing data.

**Use cases:** Quality control, read trimming, sequence alignment

## About Genomics Formats

### FASTQ
- **Content:** Raw sequencing reads + base quality scores
- **Structure:** 4 lines per read (identifier, sequence, separator, quality)
- **Quality encoding:** Phred+33 (Illumina 1.8+) or Phred+64 (older)
- **Compression:** Usually gzipped (.fastq.gz)

### SAM/BAM
- **SAM:** Human-readable text format for alignments
- **BAM:** Binary compressed version of SAM (smaller, faster)
- **Content:** Aligned reads with CIGAR strings, mapping quality
- **Tools:** samtools, Picard, GATK

### VCF/BCF
- **VCF:** Variant Call Format (text)
- **BCF:** Binary VCF (compressed)
- **Content:** SNPs, indels, structural variants
- **Fields:** Position, reference, alternate, quality, genotypes

### BED
- **Content:** Genomic regions (chromosome, start, end)
- **Uses:** Gene annotations, peaks, coverage intervals
- **Variants:** BED3, BED6, BED12

## Installation

```bash
# Python - BioPython (comprehensive bio package)
pip install biopython

# Python - pysam (SAM/BAM/VCF parsing)
pip install pysam

# Python - scikit-bio
pip install scikit-bio

# Command-line tools (bioconda recommended)
conda install -c bioconda samtools bcftools bedtools

# Or with apt (Ubuntu/Debian)
sudo apt install samtools bcftools bedtools

# FASTQ tools
conda install -c bioconda fastqc seqtk
```

## Usage Examples

### Python (BioPython - FASTQ)
```python
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Read FASTQ file
reads = list(SeqIO.parse("illumina_reads_sample.fastq", "fastq"))

print(f"Total reads: {len(reads)}")

# Examine first read
read = reads[0]
print(f"\nRead ID: {read.id}")
print(f"Sequence: {read.seq}")
print(f"Quality: {read.letter_annotations['phred_quality']}")
print(f"Length: {len(read)}")

# Calculate statistics
lengths = [len(r) for r in reads]
qualities = [np.mean(r.letter_annotations['phred_quality']) for r in reads]

print(f"\nMean read length: {np.mean(lengths):.1f} bp")
print(f"Mean quality score: {np.mean(qualities):.1f}")

# Plot quality distribution
plt.figure(figsize=(10, 6))
plt.hist(qualities, bins=20, edgecolor='black', alpha=0.7)
plt.xlabel("Mean Quality Score")
plt.ylabel("Number of Reads")
plt.title("Read Quality Distribution")
plt.grid(True, alpha=0.3)
plt.savefig("quality_distribution.png")

# Filter by quality
high_quality = [r for r in reads if np.mean(r.letter_annotations['phred_quality']) >= 30]
print(f"\nHigh quality reads (Q>=30): {len(high_quality)}")

# Write filtered reads
SeqIO.write(high_quality, "filtered_reads.fastq", "fastq")
```

### Python (pysam - SAM/BAM)
```python
import pysam
import pandas as pd

# Read SAM file
samfile = pysam.AlignmentFile("toy_alignment.sam", "r")

# Print header
print("Reference sequences:")
for ref in samfile.references:
    length = samfile.get_reference_length(ref)
    print(f"  {ref}: {length} bp")

# Iterate through alignments
alignments = []
for read in samfile:
    alignments.append({
        'name': read.query_name,
        'reference': read.reference_name,
        'position': read.reference_start,
        'mapq': read.mapping_quality,
        'flag': read.flag,
        'cigar': read.cigarstring,
        'sequence': read.query_sequence[:20] + "..." if read.query_sequence else None
    })

df = pd.DataFrame(alignments)
print(f"\nTotal alignments: {len(df)}")
print(df.head())

# Statistics
print(f"\nMapped reads: {(df['mapq'] > 0).sum()}")
print(f"Mean mapping quality: {df['mapq'].mean():.2f}")

# Count by reference
print("\nReads per reference:")
print(df['reference'].value_counts())

samfile.close()

# For BAM files (binary, must specify 'rb')
# bamfile = pysam.AlignmentFile("file.bam", "rb")
```

### Python (pysam - VCF)
```python
import pysam

# Read VCF file
vcf = pysam.VariantFile("sample_variants.vcf")

# Print header
print("VCF version:", vcf.header.version)
print("\nContigs:")
for contig in vcf.header.contigs:
    print(f"  {contig}")

print("\nSamples:")
for sample in vcf.header.samples:
    print(f"  {sample}")

# Iterate through variants
variants = []
for rec in vcf:
    variants.append({
        'chrom': rec.chrom,
        'pos': rec.pos,
        'id': rec.id,
        'ref': rec.ref,
        'alt': ','.join([str(a) for a in rec.alts]) if rec.alts else '.',
        'qual': rec.qual,
        'filter': ','.join(rec.filter.keys()) if rec.filter else 'PASS'
    })

df = pd.DataFrame(variants)
print(f"\nTotal variants: {len(df)}")
print(df.head())

# Variant type counts
df['type'] = df.apply(lambda x: 'SNP' if len(x['ref']) == 1 and len(x['alt']) == 1 else 'INDEL', axis=1)
print("\nVariant types:")
print(df['type'].value_counts())

vcf.close()
```

### Command-Line Tools

```bash
# FASTQ quality control
fastqc illumina_reads_sample.fastq

# Count FASTQ reads
echo $(($(wc -l < illumina_reads_sample.fastq) / 4))

# Convert FASTQ to FASTA
seqtk seq -a illumina_reads_sample.fastq > reads.fasta

# SAM/BAM operations
# View SAM file
samtools view toy_alignment.sam | head

# Convert SAM to BAM
samtools view -b toy_alignment.sam > alignment.bam

# Sort BAM
samtools sort alignment.bam -o alignment.sorted.bam

# Index BAM
samtools index alignment.sorted.bam

# Get alignment statistics
samtools flagstat alignment.bam
samtools stats alignment.bam

# VCF operations
# View VCF
bcftools view sample_variants.vcf

# Filter variants by quality
bcftools view -i 'QUAL>30' sample_variants.vcf

# Convert VCF to BCF (binary)
bcftools view sample_variants.vcf -O b -o variants.bcf

# Count variants
bcftools stats sample_variants.vcf

# Extract specific regions
bcftools view sample_variants.vcf chr1:1000-2000
```

## FASTQ Quality Encoding

```python
from Bio import SeqIO

def decode_quality(qual_str, encoding='phred33'):
    """Decode quality string to Phred scores"""
    offset = 33 if encoding == 'phred33' else 64
    return [ord(c) - offset for c in qual_str]

# Read and decode
for record in SeqIO.parse("illumina_reads_sample.fastq", "fastq"):
    qual = record.letter_annotations['phred_quality']
    print(f"Quality scores: {qual}")

    # Probability of error
    error_prob = [10**(-q/10) for q in qual]
    print(f"Error probability: {[f'{p:.4f}' for p in error_prob[:5]]}")
    break
```

## Common Analysis Workflows

### Read Quality Control
```python
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def analyze_fastq(filename):
    """Comprehensive FASTQ QC"""
    reads = list(SeqIO.parse(filename, "fastq"))

    # Basic stats
    print(f"Total reads: {len(reads)}")

    lengths = [len(r) for r in reads]
    mean_quals = [np.mean(r.letter_annotations['phred_quality']) for r in reads]

    print(f"Read length: {np.mean(lengths):.1f} ± {np.std(lengths):.1f} bp")
    print(f"Mean quality: {np.mean(mean_quals):.2f}")

    # Per-position quality
    max_len = max(lengths)
    pos_quals = [[] for _ in range(max_len)]

    for read in reads:
        quals = read.letter_annotations['phred_quality']
        for i, q in enumerate(quals):
            pos_quals[i].append(q)

    mean_pos_qual = [np.mean(pq) for pq in pos_quals if pq]

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(mean_pos_qual)
    ax1.axhline(y=30, color='r', linestyle='--', label='Q30')
    ax1.axhline(y=20, color='orange', linestyle='--', label='Q20')
    ax1.set_xlabel('Position in read (bp)')
    ax1.set_ylabel('Mean Quality Score')
    ax1.set_title('Quality by Position')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.hist(mean_quals, bins=30, edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Mean Read Quality')
    ax2.set_ylabel('Count')
    ax2.set_title('Read Quality Distribution')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fastq_qc.png', dpi=150)

analyze_fastq("illumina_reads_sample.fastq")
```

### Variant Filtering
```python
import pysam

def filter_variants(input_vcf, output_vcf, min_qual=30, min_depth=10):
    """Filter VCF by quality and depth"""
    vcf_in = pysam.VariantFile(input_vcf)
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    total = 0
    passed = 0

    for rec in vcf_in:
        total += 1

        # Apply filters
        if rec.qual and rec.qual >= min_qual:
            # Check depth if available
            if 'DP' in rec.info:
                if rec.info['DP'] >= min_depth:
                    vcf_out.write(rec)
                    passed += 1
            else:
                vcf_out.write(rec)
                passed += 1

    print(f"Total variants: {total}")
    print(f"Passed filters: {passed} ({100*passed/total:.1f}%)")

    vcf_in.close()
    vcf_out.close()

filter_variants("sample_variants.vcf", "filtered_variants.vcf")
```

## File Conversions

```python
from Bio import SeqIO

# FASTQ to FASTA
SeqIO.convert("illumina_reads_sample.fastq", "fastq",
              "reads.fasta", "fasta")

# Extract sequences only
with open("sequences.txt", "w") as out:
    for record in SeqIO.parse("illumina_reads_sample.fastq", "fastq"):
        out.write(f">{record.id}\n{record.seq}\n")
```

## Applications

Genomics data is used for:
- **Clinical genomics:** Disease diagnosis, personalized medicine
- **Population genetics:** Evolution, ancestry, genetic diversity
- **Cancer genomics:** Tumor profiling, mutation analysis
- **Agricultural genomics:** Crop improvement, livestock breeding
- **Metagenomics:** Microbiome analysis, environmental surveys
- **Transcriptomics:** RNA-seq, gene expression studies
- **Epigenomics:** DNA methylation, histone modifications

## Major Public Datasets

- **1000 Genomes Project:** Human genetic variation - [internationalgenome.org](http://www.internationalgenome.org/)
- **TCGA:** Cancer genomics (TB-scale) - [cancer.gov/tcga](https://www.cancer.gov/tcga)
- **gnomAD:** Human variants (>700TB) - [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/)
- **ENCODE:** Functional elements - [encodeproject.org](https://www.encodeproject.org/)
- **UK Biobank:** 500k genomes - [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/)
- **SRA:** Sequence Read Archive (PB-scale) - [ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)

## Resources

- **SAM/BAM Specification:** [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/)
- **VCF Specification:** [samtools.github.io/hts-specs/VCFv4.3.pdf](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- **BioPython Tutorial:** [biopython.org/DIST/docs/tutorial/Tutorial.html](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- **pysam Documentation:** [pysam.readthedocs.io](https://pysam.readthedocs.io/)
- **samtools:** [samtools.github.io](http://www.htslib.org/)
- **GATK Best Practices:** [gatk.broadinstitute.org](https://gatk.broadinstitute.org/)
