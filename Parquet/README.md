# Parquet Datasets

Apache Parquet is a columnar storage format optimized for analytics and machine learning workflows. While less common in traditional HPC than HDF5/NetCDF, Parquet is growing in scientific computing especially for ML applications.

## Available Datasets

*Currently building this collection. Small example datasets will be added here.*

## Shadow Datasets (Large Collections)

For Parquet datasets too large to host on GitHub, see our shadow dataset documentation:

### Computational Fluid Dynamics
- **AirfRANS** (Moderate size): [Shadow-Datasets/Computational-Fluid-Dynamics.md](../Shadow-Datasets/Computational-Fluid-Dynamics.md)
  - RANS simulations over airfoils
  - Optimized for ML pipelines
  - Pressure and velocity fields

- **DrivAerNet++** (39 TB):
  - Includes Parquet alongside HDF5/VTK
  - Automotive aerodynamics
  - 8,000 parametric car designs

## Parquet vs HDF5

| Feature | Parquet | HDF5 |
|---------|---------|------|
| **Storage** | Columnar | Row-based or chunked |
| **Use Case** | Analytics, ML | Scientific simulations |
| **Parallel I/O** | Limited | Excellent (MPI-IO) |
| **Compression** | Per-column | Per-chunk |
| **Metadata** | Schema-based | Hierarchical, flexible |
| **Ecosystem** | Spark, Pandas, Arrow | HPC tools, domain codes |

**When to use Parquet:**
- Machine learning training data
- Large-scale analytics
- Cloud-based workflows (S3, etc.)
- Selective column access

**When to use HDF5:**
- Parallel HPC applications
- Complex hierarchical data
- Scientific metadata requirements
- Multi-dimensional arrays

## Working with Parquet

### Python (pandas)
```python
import pandas as pd
import numpy as np

# Read Parquet
df = pd.read_parquet('data.parquet')

# Read specific columns only
df = pd.read_parquet('data.parquet', columns=['temperature', 'pressure'])

# Read with filters
df = pd.read_parquet('data.parquet',
                     filters=[('temperature', '>', 300)])

# Write Parquet
df = pd.DataFrame({
    'time': np.arange(1000),
    'temperature': np.random.rand(1000) * 100,
    'pressure': np.random.rand(1000) * 1000
})

df.to_parquet('output.parquet',
              compression='snappy',  # or 'gzip', 'brotli', 'zstd'
              index=False)

# Partition by column
df.to_parquet('partitioned_data/',
              partition_cols=['year', 'month'],
              compression='snappy')
```

### Python (PyArrow)
```python
import pyarrow.parquet as pq
import pyarrow as pa
import numpy as np

# Read Parquet file
table = pq.read_table('data.parquet')

# Read specific columns
table = pq.read_table('data.parquet', columns=['temperature', 'pressure'])

# Convert to pandas
df = table.to_pandas()

# Write Parquet
data = {
    'temperature': np.random.rand(1000),
    'pressure': np.random.rand(1000),
    'velocity_x': np.random.rand(1000),
    'velocity_y': np.random.rand(1000)
}

table = pa.table(data)
pq.write_table(table, 'output.parquet',
               compression='snappy',
               use_dictionary=True)

# Write with custom schema
schema = pa.schema([
    ('temperature', pa.float64()),
    ('pressure', pa.float64()),
    ('metadata', pa.string())
])

arrays = [
    pa.array(np.random.rand(1000)),
    pa.array(np.random.rand(1000)),
    pa.array(['sim1'] * 1000)
]

table = pa.Table.from_arrays(arrays, schema=schema)
pq.write_table(table, 'output.parquet')
```

### Dask (Parallel Processing)
```python
import dask.dataframe as dd

# Read large Parquet dataset
ddf = dd.read_parquet('large_dataset/*.parquet')

# Lazy operations
result = ddf[ddf['temperature'] > 300].groupby('region').mean()

# Compute result
computed = result.compute()

# Write partitioned Parquet
ddf.to_parquet('output/',
               partition_on=['year', 'month'],
               compression='snappy')
```

### PySpark
```python
from pyspark.sql import SparkSession

spark = SparkSession.builder.getOrCreate()

# Read Parquet
df = spark.read.parquet('data.parquet')

# SQL queries
df.createOrReplaceTempView('data')
result = spark.sql("""
    SELECT region, AVG(temperature) as avg_temp
    FROM data
    WHERE temperature > 300
    GROUP BY region
""")

# Write Parquet
result.write.parquet('output.parquet',
                     mode='overwrite',
                     compression='snappy')

# Partitioned write
df.write.partitionBy('year', 'month').parquet('partitioned/')
```

## Compression Options

```python
import pandas as pd

# Different compression algorithms
compressions = ['snappy', 'gzip', 'brotli', 'zstd', 'lz4']

for comp in compressions:
    df.to_parquet(f'data_{comp}.parquet', compression=comp)
```

**Compression comparison:**
- **snappy:** Fast, moderate compression (default)
- **gzip:** Slower, better compression
- **brotli:** Best compression, slower
- **zstd:** Good balance (recommended for new projects)
- **lz4:** Fastest, least compression

## Converting Between Formats

### HDF5 to Parquet
```python
import h5py
import pandas as pd

# Read from HDF5
with h5py.File('data.h5', 'r') as f:
    df = pd.DataFrame({
        'temperature': f['temperature'][:],
        'pressure': f['pressure'][:],
        'time': f['time'][:]
    })

# Write to Parquet
df.to_parquet('data.parquet', compression='snappy')
```

### Parquet to HDF5
```python
import pandas as pd

# Read Parquet
df = pd.read_parquet('data.parquet')

# Write to HDF5
df.to_hdf('data.h5', key='data', mode='w', format='table')
```

### NetCDF to Parquet (via xarray)
```python
import xarray as xr

# Read NetCDF
ds = xr.open_dataset('data.nc')

# Convert to DataFrame
df = ds.to_dataframe().reset_index()

# Write Parquet
df.to_parquet('data.parquet')
```

## Parquet in Machine Learning

### PyTorch DataLoader
```python
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader

class ParquetDataset(Dataset):
    def __init__(self, parquet_file):
        self.df = pd.read_parquet(parquet_file)

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        features = torch.tensor(row[['x', 'y', 'z']].values, dtype=torch.float32)
        label = torch.tensor(row['label'], dtype=torch.long)
        return features, label

dataset = ParquetDataset('train_data.parquet')
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
```

### TensorFlow tf.data
```python
import tensorflow as tf

def parse_parquet(filename):
    dataset = tf.data.experimental.make_batched_features_dataset(
        file_pattern=filename,
        batch_size=32,
        features={
            'temperature': tf.io.FixedLenFeature([], tf.float32),
            'pressure': tf.io.FixedLenFeature([], tf.float32),
            'label': tf.io.FixedLenFeature([], tf.int64)
        },
        reader=tf.data.experimental.make_parquet_batch_reader
    )
    return dataset

dataset = parse_parquet('train_data.parquet')
```

## Schema Management

```python
import pyarrow as pa
import pyarrow.parquet as pq

# Define schema
schema = pa.schema([
    pa.field('timestamp', pa.timestamp('ms')),
    pa.field('sensor_id', pa.int32()),
    pa.field('temperature', pa.float64(), metadata={'unit': 'celsius'}),
    pa.field('pressure', pa.float64(), metadata={'unit': 'pascal'}),
    pa.field('location', pa.struct([
        pa.field('lat', pa.float64()),
        pa.field('lon', pa.float64())
    ]))
])

# Create table with schema
data = {...}  # Your data
table = pa.Table.from_pydict(data, schema=schema)

# Write with metadata
pq.write_table(table, 'data.parquet')

# Read and inspect schema
table = pq.read_table('data.parquet')
print(table.schema)
print(table.schema.metadata)
```

## Resources

- **Apache Parquet:** [parquet.apache.org](https://parquet.apache.org/)
- **PyArrow:** [arrow.apache.org/docs/python](https://arrow.apache.org/docs/python/)
- **Pandas Parquet:** [pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_parquet.html](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_parquet.html)
- **Dask:** [docs.dask.org](https://docs.dask.org/)
