# Genomics Data Files

Comprehensive genomics datasets covering multiple formats and analysis workflows - from reference genomes and gene sequences to sequencing reads, alignments, variant calls, and expression data.

## Dataset Overview

This directory contains two complementary sets of genomics data:

### FASTA/HDF5 Datasets (8.5 MB)
Reference sequences, population genetics, and expression data in structured formats.

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| Synthetic Genome Reference | FASTA | 179 KB | 5 chromosomes, 180 Kbp total |
| Gene Sequences | FASTA | 1.1 KB | 5 annotated gene sequences |
| Variant Calls | HDF5 | 1.0 MB | 5,000 SNPs/indels × 100 samples |
| RNA-seq Expression | HDF5 | 7.3 MB | 15,000 genes × 48 samples |

### NGS Workflow Datasets (~3 KB)
Standard next-generation sequencing formats for alignment and variant analysis.

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| toy_alignment.sam | SAM | 786 bytes | Aligned sequencing reads |
| sample_variants.vcf | VCF | 2.1 KB | SNPs, insertions, deletions |
| illumina_reads_sample.fastq | FASTQ | 178 bytes | Raw reads with quality scores |

**Total:** ~8.5 MB across all datasets

---

# Part 1: FASTA/HDF5 Datasets

## 1. Synthetic Genome Reference

**File:** `synthetic_genome_reference.fasta` (179 KB)

A small reference genome with realistic GC content and chromosome structure.

**Contents:**
- 5 chromosomes (chr1-chr5)
- Total length: 180,000 base pairs
- GC content: 45-52% (varies by chromosome)
- FASTA format with 80-character line wrapping

**Chromosome Specifications:**
```
chr1: 50,000 bp, GC=45%
chr2: 40,000 bp, GC=50%
chr3: 35,000 bp, GC=48%
chr4: 30,000 bp, GC=52%
chr5: 25,000 bp, GC=47% (mitochondrial-like)
```

**Usage Example (Python):**
```python
from Bio import SeqIO

# Read FASTA file
for record in SeqIO.parse('synthetic_genome_reference.fasta', 'fasta'):
    print(f"{record.id}: {len(record.seq)} bp")
    gc_content = (record.seq.count('G') + record.seq.count('C')) / len(record.seq)
    print(f"  GC content: {gc_content:.2%}\n")
```

**Usage Example (Command Line):**
```bash
# Count sequences
grep -c "^>" synthetic_genome_reference.fasta

# Extract chromosome 1
awk '/^>chr1/,/^>/' synthetic_genome_reference.fasta | grep -v "^>chr[2-5]" > chr1.fasta

# Compute statistics with seqkit
seqkit stats synthetic_genome_reference.fasta
```

---

## 2. Gene Sequences

**File:** `gene_sequences.fasta` (1.1 KB)

Annotated gene sequences with functional descriptions.

**Contents:**
- 5 gene sequences (150-220 bp each)
- Gene annotations: ABC transporter, kinase, transcription factor, etc.
- Metadata includes organism and sequence length

**Usage Example (Python):**
```python
from Bio import SeqIO

genes = list(SeqIO.parse('gene_sequences.fasta', 'fasta'))
for gene in genes:
    print(f"Gene: {gene.id}")
    print(f"Length: {len(gene.seq)} bp")
    print(f"Description: {gene.description}\n")
```

---

## 3. Variant Calls (VCF-like HDF5)

**File:** `variant_calls.h5` (1.0 MB)

Population genetic variation data in HDF5 format.

**Contents:**
- **5,000 variants** (SNPs, insertions, deletions)
- **100 samples** with genotype calls
- Variant annotations: gene names, functional consequences
- Quality scores and allele frequencies
- Genotypes follow Hardy-Weinberg equilibrium

**HDF5 Structure:**
```
/variants/
  ├── chromosome (chr1-chr5)
  ├── position (1-based coordinates)
  ├── ref_allele
  ├── alt_allele
  ├── quality (Phred-scaled)
  ├── allele_frequency
  └── variant_type (SNP, insertion, deletion)

/genotypes/
  ├── sample_ids (SAMPLE_0001 to SAMPLE_0100)
  ├── genotypes (5000 × 100 matrix)
  │   Values: 0=ref/ref, 1=ref/alt, 2=alt/alt, -1=missing
  └── read_depth (coverage per variant per sample)

/annotations/
  ├── gene (associated gene or "intergenic")
  └── consequence (synonymous, missense, nonsense, etc.)
```

**Variant Type Distribution:**
- SNPs: 85%
- Insertions: 8%
- Deletions: 7%

**Functional Consequences:**
- Synonymous: 30%
- Missense: 25%
- Nonsense: 2%
- Intron: 20%
- Intergenic: 15%
- UTR: 8%

**Usage Example (Python):**
```python
import h5py
import numpy as np

# Read variant data
with h5py.File('variant_calls.h5', 'r') as f:
    # Get basic info
    print(f"Variants: {len(f['variants/position'][:])}")
    print(f"Samples: {f['genotypes'].attrs['n_samples']}")

    # Extract high-quality SNPs
    quality = f['variants/quality'][:]
    var_type = f['variants/variant_type'][:].astype(str)

    high_qual_snps = (quality > 30) & (var_type == 'SNP')
    print(f"High-quality SNPs: {high_qual_snps.sum()}")

    # Compute allele frequency spectrum
    af = f['variants/allele_frequency'][:]
    import matplotlib.pyplot as plt
    plt.hist(af, bins=50)
    plt.xlabel('Allele Frequency')
    plt.ylabel('Count')
    plt.title('Site Frequency Spectrum')
    plt.show()

    # Get genotypes for first 10 variants
    genotypes = f['genotypes/genotypes'][:10, :]
    sample_ids = f['genotypes/sample_ids'][:].astype(str)
```

**Usage Example (Conversion to VCF):**
```python
import h5py

with h5py.File('variant_calls.h5', 'r') as f:
    with open('variants.vcf', 'w') as vcf:
        # Write VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write variants (simplified)
        chroms = f['variants/chromosome'][:].astype(str)
        positions = f['variants/position'][:]
        refs = f['variants/ref_allele'][:].astype(str)
        alts = f['variants/alt_allele'][:].astype(str)
        quals = f['variants/quality'][:]

        for i in range(len(positions)):
            vcf.write(f"{chroms[i]}\t{positions[i]}\t.\t{refs[i]}\t{alts[i]}\t{quals[i]:.1f}\t.\t.\n")
```

---

## 4. RNA-seq Gene Expression

**File:** `rnaseq_expression.h5` (7.3 MB)

Differential gene expression data from a time-course experiment.

**Contents:**
- **15,000 genes** (GENE00001 - GENE15000)
- **48 samples** (2 conditions × 3 timepoints × 8 replicates)
- Raw read counts and TPM-normalized expression
- 1,500 differentially expressed genes between conditions

**Experimental Design:**
- **Conditions:** Control vs. Treatment
- **Timepoints:** 0, 6, 24 hours
- **Replicates:** 8 biological replicates per condition/timepoint
- **Sequencing depth:** Simulated negative binomial distribution

**HDF5 Structure:**
```
/genes/
  ├── gene_ids (GENE00001, GENE00002, ...)
  ├── gene_names
  └── length (gene length in bp)

/samples/
  ├── sample_ids (S001 - S048)
  ├── condition (control, treatment)
  ├── timepoint (0, 6, 24 hours)
  └── replicate (1-8)

/expression/
  ├── raw_counts (15000 × 48 matrix)
  └── tpm (Transcripts Per Million, 15000 × 48 matrix)
```

**Usage Example (Python):**
```python
import h5py
import numpy as np
import pandas as pd

with h5py.File('rnaseq_expression.h5', 'r') as f:
    # Create expression dataframe
    gene_ids = f['genes/gene_ids'][:].astype(str)
    sample_ids = f['samples/sample_ids'][:].astype(str)
    tpm = f['expression/tpm'][:]

    expr_df = pd.DataFrame(tpm, index=gene_ids, columns=sample_ids)

    # Get sample metadata
    conditions = f['samples/condition'][:].astype(str)
    timepoints = f['samples/timepoint'][:]

    sample_meta = pd.DataFrame({
        'sample_id': sample_ids,
        'condition': conditions,
        'timepoint': timepoints
    })

    # Find highly expressed genes
    mean_expr = expr_df.mean(axis=1)
    top_genes = mean_expr.nlargest(10)
    print("Top 10 expressed genes:")
    print(top_genes)

    # Differential expression (simple t-test)
    from scipy import stats

    control_mask = conditions == 'control'
    treatment_mask = conditions == 'treatment'

    pvalues = []
    fold_changes = []

    for i in range(len(gene_ids)):
        control_expr = tpm[i, control_mask]
        treatment_expr = tpm[i, treatment_mask]

        t_stat, pval = stats.ttest_ind(control_expr, treatment_expr)
        fc = np.log2(treatment_expr.mean() / (control_expr.mean() + 1e-6))

        pvalues.append(pval)
        fold_changes.append(fc)

    # Volcano plot
    import matplotlib.pyplot as plt
    plt.scatter(fold_changes, -np.log10(pvalues), alpha=0.3)
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot: Treatment vs Control')
    plt.axhline(-np.log10(0.05), color='r', linestyle='--', label='p=0.05')
    plt.legend()
    plt.show()
```

**Usage Example (Integration with scanpy):**
```python
import h5py
import anndata
import scanpy as sc

with h5py.File('rnaseq_expression.h5', 'r') as f:
    # Create AnnData object
    counts = f['expression/raw_counts'][:]
    gene_ids = f['genes/gene_ids'][:].astype(str)
    sample_ids = f['samples/sample_ids'][:].astype(str)

    adata = anndata.AnnData(X=counts.T)  # Transpose: samples × genes
    adata.var_names = gene_ids
    adata.obs_names = sample_ids

    # Add metadata
    adata.obs['condition'] = f['samples/condition'][:].astype(str)
    adata.obs['timepoint'] = f['samples/timepoint'][:]

# Standard scanpy workflow
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['condition', 'timepoint'])
```

---

# Part 2: NGS Workflow Datasets (FASTQ/SAM/VCF)

Standard next-generation sequencing formats used across NGS workflows, from raw reads to variant calling.

## 5. FASTQ - Raw Sequencing Reads

**File:** `illumina_reads_sample.fastq` (178 bytes)

Raw sequencing reads from Illumina platform with Phred quality scores. FASTQ is the de facto standard for storing unaligned sequencing data.

**Use cases:** Quality control, read trimming, sequence alignment

**Format Specification:**
- **Content:** Raw sequencing reads + base quality scores
- **Structure:** 4 lines per read (identifier, sequence, separator, quality)
- **Quality encoding:** Phred+33 (Illumina 1.8+) or Phred+64 (older)
- **Compression:** Usually gzipped (.fastq.gz)

**Usage Example (Python - BioPython):**
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

**Quality Control Workflow:**
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

**Command-Line Tools:**
```bash
# FASTQ quality control
fastqc illumina_reads_sample.fastq

# Count FASTQ reads
echo $(($(wc -l < illumina_reads_sample.fastq) / 4))

# Convert FASTQ to FASTA
seqtk seq -a illumina_reads_sample.fastq > reads.fasta
```

---

## 6. SAM - Sequence Alignment

**File:** `toy_alignment.sam` (786 bytes)

Example alignment file showing mapped sequencing reads. SAM is the human-readable text format for storing biological sequences aligned to a reference genome.

**Use cases:** Alignment viewing, format learning, pipeline testing

**Format Specification:**
- **SAM:** Human-readable text format for alignments
- **BAM:** Binary compressed version of SAM (smaller, faster)
- **Content:** Aligned reads with CIGAR strings, mapping quality
- **Tools:** samtools, Picard, GATK

**Usage Example (Python - pysam):**
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

**Command-Line Operations:**
```bash
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
```

---

## 7. VCF - Variant Call Format

**File:** `sample_variants.vcf` (2.1 KB)

Example variant calls showing SNPs, insertions, and deletions. VCF is the standard format for storing gene sequence variations.

**Use cases:** Variant analysis, genotype studies, population genetics

**Format Specification:**
- **VCF:** Variant Call Format (text)
- **BCF:** Binary VCF (compressed)
- **Content:** SNPs, indels, structural variants
- **Fields:** Position, reference, alternate, quality, genotypes

**Usage Example (Python - pysam):**
```python
import pysam
import pandas as pd

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

**Variant Filtering:**
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

**Command-Line Operations:**
```bash
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

---

# Installation & Tools

## Python Packages

```bash
# BioPython (comprehensive bio package)
pip install biopython

# pysam (SAM/BAM/VCF parsing)
pip install pysam

# HDF5 support
pip install h5py

# Data analysis
pip install pandas numpy scipy matplotlib

# Single-cell/expression analysis
pip install scanpy anndata

# scikit-bio
pip install scikit-bio
```

## Command-Line Tools

```bash
# Using bioconda (recommended)
conda install -c bioconda samtools bcftools bedtools

# Or with apt (Ubuntu/Debian)
sudo apt install samtools bcftools bedtools

# FASTQ tools
conda install -c bioconda fastqc seqtk

# Sequence statistics
conda install -c bioconda seqkit
```

---

# File Format Conversions

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

---

# Applications

Genomics data is used for:
- **Clinical genomics:** Disease diagnosis, personalized medicine
- **Population genetics:** Evolution, ancestry, genetic diversity
- **Cancer genomics:** Tumor profiling, mutation analysis
- **Agricultural genomics:** Crop improvement, livestock breeding
- **Metagenomics:** Microbiome analysis, environmental surveys
- **Transcriptomics:** RNA-seq, gene expression studies
- **Epigenomics:** DNA methylation, histone modifications

---

# Major Public Datasets

- **1000 Genomes Project:** Human genetic variation - [internationalgenome.org](http://www.internationalgenome.org/)
- **TCGA:** Cancer genomics (TB-scale) - [cancer.gov/tcga](https://www.cancer.gov/tcga)
- **gnomAD:** Human variants (>700TB) - [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/)
- **ENCODE:** Functional elements - [encodeproject.org](https://www.encodeproject.org/)
- **UK Biobank:** 500k genomes - [ukbiobank.ac.uk](https://www.ukbiobank.ac.uk/)
- **SRA:** Sequence Read Archive (PB-scale) - [ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)

---

# Resources

- **SAM/BAM Specification:** [samtools.github.io/hts-specs](https://samtools.github.io/hts-specs/)
- **VCF Specification:** [samtools.github.io/hts-specs/VCFv4.3.pdf](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- **BioPython Tutorial:** [biopython.org/DIST/docs/tutorial/Tutorial.html](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
- **pysam Documentation:** [pysam.readthedocs.io](https://pysam.readthedocs.io/)
- **samtools:** [samtools.github.io](http://www.htslib.org/)
- **GATK Best Practices:** [gatk.broadinstitute.org](https://gatk.broadinstitute.org/)
