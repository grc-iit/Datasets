# Genomics Datasets

This directory contains genomics datasets in FASTA and HDF5 formats for bioinformatics benchmarking and testing.

## Datasets Overview

Total size: **8.5 MB**

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| Synthetic Genome Reference | FASTA | 179 KB | 5 chromosomes, 180 Kbp total |
| Gene Sequences | FASTA | 1.1 KB | 5 annotated gene sequences |
| Variant Calls | HDF5 | 1.0 MB | 5,000 SNPs/indels × 100 samples |
| RNA-seq Expression | HDF5 | 7.3 MB | 15,000 genes × 48 samples |

---

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

**Usage Example (Integration with scanpy for single-cell-like analysis):**
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
    adata.obs_names = sample_ids
    adata.var_names = gene_ids

    # Add metadata
    adata.obs['condition'] = f['samples/condition'][:].astype(str)
    adata.obs['timepoint'] = f['samples/timepoint'][:]

    # Basic preprocessing
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # PCA and visualization
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata)
    sc.pl.pca(adata, color='condition')
```

---

## File Formats

### FASTA Format
- Standard bioinformatics sequence format
- Compatible with: BioPython, BioPerl, samtools, BLAST, etc.
- 80-character line wrapping for sequences

### HDF5 Format
- Hierarchical Data Format 5
- Efficient storage and access for large genomics data
- Compression enabled (gzip, level 6-9)
- Compatible with: h5py, PyTables, HDF5 libraries in R/MATLAB

---

## Tools and Libraries

### Python
```bash
pip install biopython h5py numpy pandas scipy
```

**Key Libraries:**
- `biopython`: FASTA parsing and sequence analysis
- `h5py`: HDF5 file access
- `pysam`: SAM/BAM/VCF file handling
- `pandas`: Data manipulation
- `scanpy`: Single-cell/bulk RNA-seq analysis

### R
```r
install.packages(c("rhdf5", "Biostrings", "DESeq2"))
```

### Command Line Tools
- `seqkit`: FASTA/FASTQ manipulation
- `samtools`: Reference genome indexing
- `bcftools`: VCF manipulation
- `h5dump`: HDF5 inspection

---

## Scientific Applications

1. **Reference Genome:**
   - Alignment algorithm testing
   - Variant calling pipeline validation
   - Genome assembly benchmarking

2. **Variant Calls:**
   - Population genetics analysis
   - GWAS (Genome-Wide Association Studies)
   - Variant effect prediction
   - Genotype imputation testing

3. **RNA-seq Expression:**
   - Differential expression analysis
   - Time-series analysis
   - Machine learning for gene expression
   - Clustering and classification

---

## Citation

```
GRC Datasets Repository - Genomics Collection
Created: October 2024
Formats: FASTA, HDF5
Repository: https://github.com/grc-iit/Datasets
```

---

## License

These datasets are provided for educational, research, and benchmarking purposes.
