# Seismology Datasets

This directory contains seismological datasets in HDF5 format simulating earthquake waveforms and seismic analysis.

## Datasets Overview

Total size: **35 MB**

| Dataset | Format | Size | Description |
|---------|--------|------|-------------|
| Earthquake Waveforms | HDF5 | 9.2 MB | 12 stations × 3 components (Z,N,E) |
| Ambient Noise | HDF5 | 26 MB | 24 hours continuous seismic noise |
| Receiver Functions | HDF5 | 449 KB | 50 teleseismic RFs + stack |

---

## 1. Earthquake Waveforms

**File:** `earthquake_waveforms.h5` (9.2 MB)

Three-component seismograms from a magnitude 5.8 earthquake recorded at 12 stations.

**Event Information:**
- **Event ID:** EQ20240815_123045
- **Origin time:** 2024-08-15T12:30:45Z
- **Location:** Southern California (35.5°N, -118.2°W)
- **Depth:** 12.5 km
- **Magnitude:** 5.8 Mw

**Recording Parameters:**
- **Stations:** 12 (network XX)
- **Components:** Z (vertical), N (north), E (east)
- **Sample rate:** 100 Hz
- **Duration:** 180 seconds (3 minutes)
- **Distance range:** ~30-150 km from epicenter

**Seismic Phases:**
- **P-wave arrivals:** 5-25 seconds (depends on distance)
- **S-wave arrivals:** 8.5-43 seconds
- **Surface waves:** Later arrivals with larger amplitude
- **Background noise:** ~50 counts RMS

**HDF5 Structure:**
```
/event/
  attributes: event_id, origin_time, latitude, longitude, depth_km, magnitude

/stations/
  /STA01/ through /STA12/
    attributes: network, station, lat, lon, elevation_m, distance_km, azimuth
    ├── Z/ (vertical component)
    │   ├── data (18000 samples)
    │   ├── time (seconds)
    │   └── attributes: sample_rate, p_arrival_sec, s_arrival_sec
    ├── N/ (north component)
    └── E/ (east component)
```

**Usage Example (Python):**
```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Read earthquake data
with h5py.File('earthquake_waveforms.h5', 'r') as f:
    # Event info
    event = f['event']
    print(f"Event: {event.attrs['event_id']}")
    print(f"Magnitude: {event.attrs['magnitude']}")
    print(f"Location: {event.attrs['location']}")

    # Plot three-component seismogram for station 1
    sta = f['stations/STA01']
    print(f"\nStation: {sta.attrs['station']}")
    print(f"Distance: {sta.attrs['distance_km']:.1f} km")

    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)

    for ax, comp in zip(axes, ['Z', 'N', 'E']):
        time = sta[comp]['time'][:]
        data = sta[comp]['data'][:]

        ax.plot(time, data, 'k-', linewidth=0.5)
        ax.set_ylabel(f'{comp} (counts)')
        ax.grid(alpha=0.3)

        # Mark phase arrivals
        p_arrival = sta[comp].attrs['p_arrival_sec']
        s_arrival = sta[comp].attrs['s_arrival_sec']

        ax.axvline(p_arrival, color='b', linestyle='--', alpha=0.7, label='P')
        ax.axvline(s_arrival, color='r', linestyle='--', alpha=0.7, label='S')

        if comp == 'Z':
            ax.legend()

    axes[-1].set_xlabel('Time (seconds)')
    plt.suptitle(f"Three-Component Seismogram - {sta.attrs['station']}")
    plt.tight_layout()
    plt.show()
```

**Travel Time Analysis:**
```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File('earthquake_waveforms.h5', 'r') as f:
    stations = f['stations']

    distances = []
    p_times = []
    s_times = []

    for sta_name in stations.keys():
        sta = stations[sta_name]
        distances.append(sta.attrs['distance_km'])
        p_times.append(sta['Z'].attrs['p_arrival_sec'])
        s_times.append(sta['Z'].attrs['s_arrival_sec'])

    distances = np.array(distances)
    p_times = np.array(p_times)
    s_times = np.array(s_times)

    # Plot travel time curves
    plt.figure(figsize=(10, 6))
    plt.plot(distances, p_times, 'bo', label='P-wave', markersize=8)
    plt.plot(distances, s_times, 'ro', label='S-wave', markersize=8)

    # Fit velocity models
    p_velocity = np.polyfit(distances, p_times, 1)[0]**-1
    s_velocity = np.polyfit(distances, s_times, 1)[0]**-1

    plt.plot(distances, distances/6.0, 'b--', alpha=0.5, label=f'P: 6.0 km/s')
    plt.plot(distances, distances/3.5, 'r--', alpha=0.5, label=f'S: 3.5 km/s')

    plt.xlabel('Distance (km)')
    plt.ylabel('Travel Time (s)')
    plt.title('Seismic Phase Travel Times')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

    print(f"Measured P-wave velocity: {p_velocity:.2f} km/s")
    print(f"Measured S-wave velocity: {s_velocity:.2f} km/s")
```

**Spectral Analysis:**
```python
from scipy import signal

with h5py.File('earthquake_waveforms.h5', 'r') as f:
    sta = f['stations/STA01']
    data = sta['Z']['data'][:]
    fs = sta['Z'].attrs['sample_rate']

    # Compute spectrogram
    f_spec, t_spec, Sxx = signal.spectrogram(data, fs, nperseg=512)

    plt.figure(figsize=(12, 6))
    plt.pcolormesh(t_spec, f_spec, 10*np.log10(Sxx), shading='gouraud',
                   cmap='viridis', vmax=60)
    plt.ylabel('Frequency (Hz)')
    plt.xlabel('Time (s)')
    plt.title('Seismogram Spectrogram')
    plt.colorbar(label='Power (dB)')
    plt.ylim(0, 20)
    plt.show()
```

---

## 2. Ambient Seismic Noise

**File:** `ambient_noise_continuous.h5` (26 MB)

Continuous seismic background noise for 24 hours, typical of ambient noise cross-correlation studies.

**Recording Parameters:**
- **Station:** BKL (Berkeley Lab)
- **Location:** 37.8744°N, -122.2597°W, 120 m elevation
- **Sample rate:** 40 Hz
- **Duration:** 24 hours (86,400 seconds)
- **Total samples:** 3,456,000

**Noise Characteristics:**
- **Primary microseism:** ~0.07 Hz (ocean wave interactions)
- **Secondary microseism:** ~0.15 Hz (dominant peak)
- **Cultural noise:** 1-10 Hz (traffic, industry)
- **Diurnal variation:** Cultural noise stronger during daytime

**HDF5 Structure:**
```
/station_BKL/
  attributes: network, station, latitude, longitude, elevation_m
  ├── hour_00/ through hour_23/
  │   ├── data (144000 samples per hour)
  │   └── attributes: start_time, sample_rate, duration_sec
```

**Usage Example (Python):**
```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Read one hour of data
with h5py.File('ambient_noise_continuous.h5', 'r') as f:
    station = f['station_BKL']

    hour_data = station['hour_12']['data'][:]
    fs = station['hour_12'].attrs['sample_rate']

    # Plot time series (first 10 minutes)
    time = np.arange(len(hour_data[:24000])) / fs / 60  # minutes

    plt.figure(figsize=(14, 4))
    plt.plot(time, hour_data[:24000], 'k-', linewidth=0.3)
    plt.xlabel('Time (minutes)')
    plt.ylabel('Amplitude (counts)')
    plt.title('Ambient Seismic Noise - Hour 12 (Noon)')
    plt.grid(alpha=0.3)
    plt.show()
```

**Power Spectral Density:**
```python
from scipy import signal

with h5py.File('ambient_noise_continuous.h5', 'r') as f:
    # Concatenate several hours
    data_list = []
    for hour in range(0, 24, 6):  # Every 6 hours
        hour_data = f[f'station_BKL/hour_{hour:02d}/data'][:]
        data_list.append(hour_data)

    fs = f['station_BKL/hour_00'].attrs['sample_rate']

    # Compute PSD for each time window
    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (hour_idx, data) in enumerate(zip([0, 6, 12, 18], data_list)):
        f_psd, psd = signal.welch(data, fs, nperseg=4096)

        ax.loglog(f_psd, np.sqrt(psd), label=f'Hour {hour_idx:02d}')

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (counts/√Hz)')
    ax.set_title('Ambient Noise Power Spectral Density')
    ax.legend()
    ax.grid(which='both', alpha=0.3)

    # Mark microseism peaks
    ax.axvline(0.07, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0.15, color='gray', linestyle='--', alpha=0.5)
    ax.text(0.07, ax.get_ylim()[1]*0.5, 'Primary\nmicroseism', ha='center')
    ax.text(0.15, ax.get_ylim()[1]*0.5, 'Secondary\nmicroseism', ha='center')

    plt.show()
```

**Cross-Correlation (Ambient Noise Seismology):**
```python
# Cross-correlate different hours to extract Green's function
with h5py.File('ambient_noise_continuous.h5', 'r') as f:
    data1 = f['station_BKL/hour_00/data'][:]
    data2 = f['station_BKL/hour_12/data'][:]

    # Normalize
    data1 = data1 / np.std(data1)
    data2 = data2 / np.std(data2)

    # Cross-correlation
    correlation = np.correlate(data1, data2, mode='full')
    lags = np.arange(-len(data1)+1, len(data1))
    time_lags = lags / fs

    # Plot central part
    center = len(correlation) // 2
    window = 4000

    plt.figure(figsize=(12, 5))
    plt.plot(time_lags[center-window:center+window],
             correlation[center-window:center+window])
    plt.xlabel('Time lag (seconds)')
    plt.ylabel('Cross-correlation')
    plt.title('Ambient Noise Cross-Correlation')
    plt.grid(alpha=0.3)
    plt.show()
```

---

## 3. Receiver Functions

**File:** `receiver_functions.h5` (449 KB)

Teleseismic receiver functions for crustal structure analysis.

**Method:**
- **Technique:** P-to-S converted phases from teleseismic earthquakes
- **Purpose:** Imaging crustal discontinuities (esp. Moho)
- **Station:** ABC (40.5°N, -105.2°W, 1650 m)
- **Events:** 50 teleseismic earthquakes (30-90° distance)
- **Frequency:** 10 Hz sampling

**Expected Arrivals:**
- **Direct P:** t = 0 s (reference)
- **Ps (Moho conversion):** ~3.5 s (for 35 km depth)
- **PpPs (Moho multiple):** ~7 s
- **Mid-crustal discontinuity:** ~2 s

**HDF5 Structure:**
```
/station_ABC/
  attributes: network, station, latitude, longitude, elevation_m
  ├── event_001/ through event_050/
  │   ├── radial_rf (receiver function)
  │   ├── time (-5 to +25 seconds)
  │   └── attributes: distance_deg, back_azimuth, magnitude
  └── stack/
      ├── radial_rf_stack (stacked RF)
      ├── time
      └── attributes: n_events_stacked
```

**Usage Example (Python):**
```python
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Read and plot individual RFs
with h5py.File('receiver_functions.h5', 'r') as f:
    station = f['station_ABC']

    # Plot first 10 RFs
    fig, ax = plt.subplots(figsize=(10, 8))

    for i in range(10):
        event = station[f'event_{i:03d}']
        time = event['time'][:]
        rf = event['radial_rf'][:]

        # Offset for display
        ax.plot(time, rf + i*0.3, 'k-', linewidth=0.8)

    ax.set_xlabel('Time (seconds)')
    ax.set_ylabel('Receiver Function (offset for display)')
    ax.set_title('Individual Receiver Functions')
    ax.axvline(0, color='b', linestyle='--', alpha=0.3, label='P arrival')
    ax.axvline(3.5, color='r', linestyle='--', alpha=0.3, label='Ps (Moho)')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.show()
```

**Stacked Receiver Function:**
```python
with h5py.File('receiver_functions.h5', 'r') as f:
    stack = f['station_ABC/stack']
    time = stack['time'][:]
    rf_stack = stack['radial_rf_stack'][:]

    n_events = stack.attrs['n_events_stacked']

    plt.figure(figsize=(12, 6))
    plt.plot(time, rf_stack, 'k-', linewidth=2)
    plt.fill_between(time, 0, rf_stack, where=(rf_stack > 0),
                    color='red', alpha=0.3, label='Positive (Ps)')
    plt.fill_between(time, 0, rf_stack, where=(rf_stack < 0),
                    color='blue', alpha=0.3, label='Negative')

    plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(3.5, color='red', linestyle='--', alpha=0.7, label='Moho Ps')
    plt.axvline(7.0, color='orange', linestyle='--', alpha=0.7, label='PpPs')

    plt.xlabel('Time (seconds)')
    plt.ylabel('Amplitude')
    plt.title(f'Stacked Receiver Function (N={n_events})')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xlim(-2, 15)
    plt.show()
```

**H-κ Stacking (Moho Depth and Vp/Vs):**
```python
# Grid search over crustal thickness (H) and Vp/Vs ratio (κ)
with h5py.File('receiver_functions.h5', 'r') as f:
    station = f['station_ABC']

    # Parameter ranges
    H_range = np.linspace(25, 45, 50)  # km
    kappa_range = np.linspace(1.65, 1.85, 50)  # Vp/Vs

    stack_grid = np.zeros((len(H_range), len(kappa_range)))

    # Simplified H-κ stacking
    for event_idx in range(50):
        event = station[f'event_{event_idx:03d}']
        time = event['time'][:]
        rf = event['radial_rf'][:]

        # For each H, κ pair, predict arrival times and stack
        for i, H in enumerate(H_range):
            for j, kappa in enumerate(kappa_range):
                # Predict Ps time (simplified)
                Vp = 6.3  # km/s
                Vs = Vp / kappa
                Ps_time = H * (np.sqrt(1/Vs**2 - 0.06**2) - np.sqrt(1/Vp**2 - 0.06**2))

                # Find RF amplitude at predicted time
                idx = np.argmin(np.abs(time - Ps_time))
                if 0 <= idx < len(rf):
                    stack_grid[i, j] += rf[idx]

    # Plot H-κ stack
    plt.figure(figsize=(10, 8))
    plt.contourf(kappa_range, H_range, stack_grid, levels=20, cmap='RdBu_r')
    plt.colorbar(label='Stacked Amplitude')
    plt.xlabel('Vp/Vs (κ)')
    plt.ylabel('Crustal Thickness (km)')
    plt.title('H-κ Stacking Analysis')

    # Find maximum
    max_idx = np.unravel_index(stack_grid.argmax(), stack_grid.shape)
    best_H = H_range[max_idx[0]]
    best_kappa = kappa_range[max_idx[1]]
    plt.plot(best_kappa, best_H, 'k*', markersize=20, label=f'Best: H={best_H:.1f} km, κ={best_kappa:.2f}')
    plt.legend()
    plt.show()

    print(f"Estimated crustal thickness: {best_H:.1f} km")
    print(f"Estimated Vp/Vs ratio: {best_kappa:.2f}")
```

---

## File Format

All datasets use **HDF5** (Hierarchical Data Format 5) with:
- **Compression:** gzip (levels 6-9)
- **Chunking:** Enabled for time-series data
- **Metadata:** Extensive attributes at all levels

---

## Tools and Software

### Python
```bash
pip install h5py numpy scipy matplotlib obspy
```

**Key Libraries:**
- `h5py`: HDF5 file I/O
- `obspy`: Seismological data processing
- `scipy.signal`: Signal processing
- `numpy`: Numerical operations

### ObsPy Integration
```python
from obspy import Stream, Trace
import h5py
import numpy as np

# Convert HDF5 to ObsPy Stream
with h5py.File('earthquake_waveforms.h5', 'r') as f:
    sta = f['stations/STA01']

    stream = Stream()
    for comp in ['Z', 'N', 'E']:
        data = sta[comp]['data'][:]
        sr = sta[comp].attrs['sample_rate']

        tr = Trace(data=data)
        tr.stats.sampling_rate = sr
        tr.stats.station = sta.attrs['station']
        tr.stats.channel = f'HH{comp}'

        stream += tr

# Now use ObsPy functions
stream.filter('bandpass', freqmin=1.0, freqmax=10.0)
stream.plot()
```

### Command Line
```bash
# Inspect HDF5 structure
h5dump -H earthquake_waveforms.h5

# List all datasets
h5ls -r earthquake_waveforms.h5
```

---

## Scientific Applications

1. **Earthquake Waveforms:**
   - Seismic phase identification
   - Travel time analysis
   - Earthquake location algorithms
   - Ground motion modeling

2. **Ambient Noise:**
   - Ambient noise cross-correlation
   - Surface wave tomography
   - Noise source location
   - Long-term monitoring

3. **Receiver Functions:**
   - Crustal thickness estimation
   - Vp/Vs ratio determination
   - Subsurface imaging
   - H-κ stacking analysis

---

## References

### Standards
- **SEED/miniSEED:** https://www.fdsn.org/seed_manual/
- **SAC Format:** https://ds.iris.edu/files/sac-manual/

### Data Centers
- **IRIS DMC:** https://ds.iris.edu/
- **FDSN:** https://www.fdsn.org/
- **ISC:** http://www.isc.ac.uk/

### Software
- **ObsPy:** https://docs.obspy.org/
- **SAC:** https://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/

---

## Citation

```
GRC Datasets Repository - Seismology Collection
Created: October 2024
Format: HDF5
Repository: https://github.com/grc-iit/Datasets
```

---

## License

These datasets are provided for educational, research, and benchmarking purposes.
