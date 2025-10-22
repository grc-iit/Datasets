# NYC Taxi Trip Data - January 2024

**Format:** Apache Parquet
**Size:** 48 MB
**Source:** NYC Taxi & Limousine Commission
**Access:** [nyc.gov/tlc](https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page)

## File

### nyc_taxi_jan2024_sample.parquet
Yellow taxi trip records for January 2024 from New York City's official TLC database.

**Records:** ~2.9 million trips
**Coverage:** All yellow taxi trips in Manhattan and boroughs
**Time Period:** January 1-31, 2024

## About NYC Taxi Data

The NYC Taxi & Limousine Commission (TLC) publishes detailed trip-level data for yellow and green taxis, for-hire vehicles (Uber/Lyft), and high-volume for-hire services. This open dataset has become a standard benchmark for:
- Big data processing systems
- Geospatial analytics
- Time-series analysis
- Data engineering tutorials

## Data Schema

```python
import pandas as pd

df = pd.read_parquet('nyc_taxi_jan2024_sample.parquet')
print(df.dtypes)
```

**Columns:**
- `VendorID`: Provider identifier (1=Creative Mobile, 2=VeriFone)
- `tpep_pickup_datetime`: Trip start timestamp
- `tpep_dropoff_datetime`: Trip end timestamp
- `passenger_count`: Number of passengers
- `trip_distance`: Distance in miles
- `RatecodeID`: Rate code (1=Standard, 2=JFK, 3=Newark, etc.)
- `store_and_fwd_flag`: Whether trip was stored before sending to vendor
- `PULocationID`: Pickup taxi zone ID
- `DOLocationID`: Dropoff taxi zone ID
- `payment_type`: Payment method (1=Credit card, 2=Cash, etc.)
- `fare_amount`: Base fare
- `extra`: Miscellaneous extras
- `mta_tax`: MTA tax ($0.50)
- `tip_amount`: Tip amount (credit cards only)
- `tolls_amount`: Total tolls
- `improvement_surcharge`: Improvement surcharge ($0.30)
- `total_amount`: Total charge to passengers
- `congestion_surcharge`: Congestion surcharge ($2.50 below 96th St)
- `Airport_fee`: Airport fee ($1.25 for LGA/JFK pickups)

## Usage Examples

### Python (pandas)
```python
import pandas as pd
import matplotlib.pyplot as plt

# Read data
df = pd.read_parquet('nyc_taxi_jan2024_sample.parquet')

print(f"Total trips: {len(df):,}")
print(f"Date range: {df['tpep_pickup_datetime'].min()} to {df['tpep_pickup_datetime'].max()}")

# Basic statistics
print("\nTrip Statistics:")
print(df[['trip_distance', 'fare_amount', 'tip_amount', 'total_amount']].describe())

# Average tip percentage (credit cards only)
credit_card_tips = df[df['payment_type'] == 1].copy()
credit_card_tips['tip_pct'] = credit_card_tips['tip_amount'] / credit_card_tips['fare_amount'] * 100
print(f"\nAverage tip percentage: {credit_card_tips['tip_pct'].mean():.1f}%")

# Trips per day
df['date'] = pd.to_datetime(df['tpep_pickup_datetime']).dt.date
daily_trips = df.groupby('date').size()
daily_trips.plot(figsize=(12, 4), title='Daily Taxi Trips - January 2024')
plt.ylabel('Number of Trips')
plt.show()

# Peak hours
df['hour'] = pd.to_datetime(df['tpep_pickup_datetime']).dt.hour
hourly_trips = df.groupby('hour').size()
hourly_trips.plot(kind='bar', figsize=(10, 5), title='Trips by Hour of Day')
plt.xlabel('Hour of Day')
plt.ylabel('Number of Trips')
plt.show()
```

### Python (DuckDB - SQL on Parquet)
```python
import duckdb

# Query Parquet directly with SQL
con = duckdb.connect()

# Top 10 most popular routes
query = """
SELECT
    PULocationID,
    DOLocationID,
    COUNT(*) as num_trips,
    AVG(trip_distance) as avg_distance,
    AVG(total_amount) as avg_fare
FROM 'nyc_taxi_jan2024_sample.parquet'
GROUP BY PULocationID, DOLocationID
ORDER BY num_trips DESC
LIMIT 10
"""

result = con.execute(query).df()
print(result)

# Average fare by hour and payment type
query = """
SELECT
    HOUR(tpep_pickup_datetime) as hour,
    payment_type,
    COUNT(*) as trips,
    AVG(total_amount) as avg_fare
FROM 'nyc_taxi_jan2024_sample.parquet'
GROUP BY hour, payment_type
ORDER BY hour, payment_type
"""

hourly_by_payment = con.execute(query).df()
```

### Python (Polars - Fast DataFrame Library)
```python
import polars as pl

# Read with Polars (faster than pandas)
df = pl.read_parquet('nyc_taxi_jan2024_sample.parquet')

# Calculate trip duration
df = df.with_columns([
    (pl.col('tpep_dropoff_datetime') - pl.col('tpep_pickup_datetime'))
    .dt.total_seconds().alias('duration_seconds')
])

# Group by pickup hour
hourly_stats = (
    df.group_by(pl.col('tpep_pickup_datetime').dt.hour().alias('hour'))
    .agg([
        pl.count().alias('num_trips'),
        pl.mean('trip_distance').alias('avg_distance'),
        pl.mean('fare_amount').alias('avg_fare'),
        pl.mean('duration_seconds').alias('avg_duration_sec')
    ])
    .sort('hour')
)

print(hourly_stats)
```

### PySpark (Big Data Processing)
```python
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, hour, avg, count

spark = SparkSession.builder.appName("TaxiAnalysis").getOrCreate()

# Read Parquet
df = spark.read.parquet('nyc_taxi_jan2024_sample.parquet')

# Analysis
hourly_stats = (
    df.withColumn('hour', hour('tpep_pickup_datetime'))
    .groupBy('hour')
    .agg(
        count('*').alias('num_trips'),
        avg('trip_distance').alias('avg_distance'),
        avg('total_amount').alias('avg_fare')
    )
    .orderBy('hour')
)

hourly_stats.show()
```

## Analysis Ideas

1. **Temporal patterns:** Rush hour analysis, weekend vs weekday
2. **Geospatial:** Popular routes, taxi zones heatmaps
3. **Revenue analysis:** Fare distribution, tipping behavior
4. **Demand forecasting:** Predict trip volumes
5. **Anomaly detection:** Unusual trips or pricing
6. **Payment patterns:** Cash vs credit card trends

## Benchmark Performance

This dataset is commonly used for benchmarking:

```python
import time
import pandas as pd
import polars as pl
import duckdb

# Pandas
start = time.time()
df_pandas = pd.read_parquet('nyc_taxi_jan2024_sample.parquet')
result = df_pandas.groupby('PULocationID')['total_amount'].mean()
pandas_time = time.time() - start

# Polars
start = time.time()
df_polars = pl.read_parquet('nyc_taxi_jan2024_sample.parquet')
result = df_polars.group_by('PULocationID').agg(pl.mean('total_amount'))
polars_time = time.time() - start

# DuckDB
start = time.time()
result = duckdb.query("""
    SELECT PULocationID, AVG(total_amount)
    FROM 'nyc_taxi_jan2024_sample.parquet'
    GROUP BY PULocationID
""").df()
duckdb_time = time.time() - start

print(f"Pandas: {pandas_time:.2f}s")
print(f"Polars: {polars_time:.2f}s")
print(f"DuckDB: {duckdb_time:.2f}s")
```

## Taxi Zone Mapping

Download the taxi zone shapefile for geospatial analysis:
```bash
wget https://d37ci6vzurychx.cloudfront.net/misc/taxi_zones.zip
```

## Related Datasets

Available from the same source:
- **Green taxi trips:** Outer boroughs
- **For-hire vehicle trips:** Uber, Lyft, etc.
- **Historical data:** 2009-present
- **All months:** Full time series

## Data Dictionary

Full data dictionary: [nyc.gov/assets/tlc/downloads/pdf/data_dictionary_trip_records_yellow.pdf](https://www.nyc.gov/assets/tlc/downloads/pdf/data_dictionary_trip_records_yellow.pdf)

## References

- **NYC TLC:** [nyc.gov/tlc](https://www.nyc.gov/site/tlc/about/tlc-trip-record-data.page)
- **Download all months:** [d37ci6vzurychx.cloudfront.net/trip-data](https://d37ci6vzurychx.cloudfront.net/trip-data/)
