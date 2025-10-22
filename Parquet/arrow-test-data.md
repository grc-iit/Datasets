# Apache Arrow Test Data - All Types Sample

**Format:** Apache Parquet
**Size:** 287 KB
**Source:** Apache Arrow Project
**Access:** [github.com/apache/arrow-testing](https://github.com/apache/arrow-testing)

## File

### arrow_alltypes_sample.parquet
Sample Parquet file demonstrating all supported data types in the Parquet format. Used for testing and validation of Parquet readers/writers across different languages and libraries.

## About Apache Arrow Testing Data

The Apache Arrow project maintains a collection of test files to ensure compatibility across implementations in C++, Python, R, Java, JavaScript, and other languages. These files serve as:
- **Format validation:** Ensure correct Parquet specification implementation
- **Type coverage:** Test all supported data types
- **Interoperability:** Verify cross-language compatibility
- **Benchmarking:** Standard test cases for performance testing

## Data Schema

```python
import pyarrow.parquet as pq

table = pq.read_table('arrow_alltypes_sample.parquet')
print(table.schema)
```

**Data Types Included:**
- **Integers:** int8, int16, int32, int64, uint8, uint16, uint32, uint64
- **Floating point:** float32, float64
- **Boolean:** bool
- **String:** string (UTF-8)
- **Binary:** binary
- **Timestamp:** timestamp with various resolutions
- **Date:** date32, date64
- **Time:** time32, time64
- **Decimal:** decimal128
- **Lists:** list types
- **Structs:** nested structures
- **Null:** null type

## Usage Examples

### Python (PyArrow)
```python
import pyarrow.parquet as pq
import pyarrow as pa

# Read the file
table = pq.read_table('arrow_alltypes_sample.parquet')

print(f"Number of rows: {len(table)}")
print(f"Number of columns: {len(table.schema)}")
print(f"\nColumn names: {table.column_names}")
print(f"\nSchema:")
print(table.schema)

# Convert to pandas
df = table.to_pandas()
print(f"\nDataFrame shape: {df.shape}")
print(df.head())

# Access specific columns
if 'id' in table.column_names:
    ids = table.column('id').to_pylist()
    print(f"\nFirst 5 IDs: {ids[:5]}")

# Check data types
for i, field in enumerate(table.schema):
    print(f"{field.name}: {field.type}")
```

### Python (pandas)
```python
import pandas as pd

# Read as DataFrame
df = pd.read_parquet('arrow_alltypes_sample.parquet')

# Inspect dtypes
print("Data types:")
print(df.dtypes)

# Basic statistics
print("\nNumeric columns statistics:")
print(df.describe())

# Check for missing values
print("\nMissing values:")
print(df.isnull().sum())
```

### Python (Polars)
```python
import polars as pl

# Read with Polars
df = pl.read_parquet('arrow_alltypes_sample.parquet')

print(df.schema)
print(df.head())

# Get data type information
for col in df.columns:
    print(f"{col}: {df[col].dtype}")
```

### Inspect Parquet Metadata
```python
import pyarrow.parquet as pq

# Read metadata
parquet_file = pq.ParquetFile('arrow_alltypes_sample.parquet')

print("Metadata:")
print(f"Created by: {parquet_file.metadata.created_by}")
print(f"Number of row groups: {parquet_file.metadata.num_row_groups}")
print(f"Number of rows: {parquet_file.metadata.num_rows}")
print(f"Number of columns: {parquet_file.metadata.num_columns}")

# Column statistics
for i in range(parquet_file.metadata.num_row_groups):
    rg = parquet_file.metadata.row_group(i)
    print(f"\nRow Group {i}:")
    for col in range(rg.num_columns):
        col_meta = rg.column(col)
        print(f"  {col_meta.path_in_schema}: "
              f"size={col_meta.total_compressed_size} bytes, "
              f"encoding={col_meta.encodings}")
```

### Command Line (parquet-tools)
```bash
# Install parquet-tools
pip install parquet-tools

# Show schema
parquet-tools schema arrow_alltypes_sample.parquet

# Show sample data
parquet-tools head arrow_alltypes_sample.parquet

# Show metadata
parquet-tools meta arrow_alltypes_sample.parquet

# Convert to JSON
parquet-tools cat arrow_alltypes_sample.parquet

# Count rows
parquet-tools rowcount arrow_alltypes_sample.parquet
```

## Testing Parquet Compatibility

This file is useful for testing Parquet implementations:

```python
import pyarrow.parquet as pq
import pyarrow as pa

# Read the test file
original = pq.read_table('arrow_alltypes_sample.parquet')

# Write it back with different options
pq.write_table(
    original,
    'test_output.parquet',
    compression='snappy',      # or 'gzip', 'brotli', 'zstd', 'lz4'
    use_dictionary=True,
    write_statistics=True,
    version='2.6'              # Parquet format version
)

# Read back and verify
rewritten = pq.read_table('test_output.parquet')

# Check if data is identical
assert original.equals(rewritten), "Data mismatch after round-trip!"
print("Round-trip test passed!")

# Compare schemas
print("Original schema:", original.schema)
print("Rewritten schema:", rewritten.schema)
```

## Type Conversion Examples

```python
import pyarrow.parquet as pq
import pyarrow as pa
import pandas as pd

# Read table
table = pq.read_table('arrow_alltypes_sample.parquet')

# Convert to different formats
pandas_df = table.to_pandas()
python_dict = table.to_pydict()
numpy_dict = {col: table.column(col).to_numpy() for col in table.column_names}

# Convert back to Arrow
table_from_pandas = pa.Table.from_pandas(pandas_df)
table_from_dict = pa.Table.from_pydict(python_dict)

# Check roundtrip
assert table.schema.equals(table_from_pandas.schema)
```

## Compression Comparison

```python
import pyarrow.parquet as pq
import os

table = pq.read_table('arrow_alltypes_sample.parquet')

compressions = ['none', 'snappy', 'gzip', 'brotli', 'zstd', 'lz4']

for comp in compressions:
    filename = f'test_{comp}.parquet'
    pq.write_table(table, filename, compression=comp)
    size = os.path.getsize(filename)
    print(f"{comp:10s}: {size:8,} bytes")
    os.remove(filename)
```

## Applications

This test file is valuable for:

1. **Library development:** Testing Parquet readers/writers
2. **Format validation:** Ensuring spec compliance
3. **Performance benchmarking:** Comparing different implementations
4. **Education:** Learning Parquet data types and structure
5. **Interoperability testing:** Cross-language data exchange
6. **Compression testing:** Evaluating compression algorithms

## Schema Evolution Testing

```python
import pyarrow.parquet as pq
import pyarrow as pa

# Read original
original = pq.read_table('arrow_alltypes_sample.parquet')

# Add a new column
extended = original.append_column(
    'new_column',
    pa.array([i * 10 for i in range(len(original))])
)

# Write with new schema
pq.write_table(extended, 'extended.parquet')

# Read back (demonstrates schema evolution)
reread = pq.read_table('extended.parquet')
print("New schema:", reread.schema)
```

## References

- **Apache Arrow:** [arrow.apache.org](https://arrow.apache.org/)
- **Arrow Testing Repo:** [github.com/apache/arrow-testing](https://github.com/apache/arrow-testing)
- **Parquet Format:** [parquet.apache.org/docs](https://parquet.apache.org/docs/)
- **PyArrow Docs:** [arrow.apache.org/docs/python](https://arrow.apache.org/docs/python/)
