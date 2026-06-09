# Performance of Wireless Digital Communication System Corrupted by Noise

**Course:** ECE 4104 - Wireless Engineering Laboratory  
**Implementation Language:** MATLAB

## Project Overview

This project simulates a digital communication system utilizing **FSK (Frequency Shift Keying)** modulation to analyze system performance under AWGN (Additive White Gaussian Noise) channel conditions. The simulation evaluates key performance metrics including Bit Error Rate (BER), Packet Error Rate (PER), and compares measured results against theoretical predictions.

## System Architecture

### Modulation Scheme
- **Type:** FSK (Frequency Shift Keying) with 2 frequency tones
- **Frequency Separation:** fs/4 (where fs is sampling frequency)
- **Amplitude Control:** Variable (2V to 5V)

### Channel Model
- **Noise Type:** AWGN (Additive White Gaussian Noise)
- **SNR Range:** 0 to 12 dB

## Simulation Parameters

### Variable Parameters Tested
| Parameter | Values | Unit |
|-----------|--------|------|
| Bit Rate | 2, 3, 4, 5 | bits/sec |
| Signal Amplitude | 2, 3, 4, 5 | Volts |
| Sampling Frequency | 12, 16, 20, 24 | Hz |
| SNR (Signal-to-Noise Ratio) | 0, 2, 4, 6, 8, 10, 12 | dB |

### Fixed Test Parameters
- **Packets per Test:** 10 packets
- **Bits per Packet:** 10 bits
- **Total Bits per Simulation:** 100 bits
- **Data Generation:** Pseudorandom binary sequence with configurable seeds

## Performance Metrics

1. **Bit Error Rate (BER):** Measured and theoretical comparison
   - Theoretical BER for FSK: BER = Q(√(2·SNR))
   - Q-function based threshold detection

2. **Packet Error Rate (PER):** Packet-level error measurement
   - Packets with one or more bit errors are counted as erroneous

3. **Signal-to-Noise Ratio (SNR):** Calculated from signal amplitude and noise amplitude

## Simulation Analysis

The project generates four comprehensive analysis plots:

1. **BER vs Bit Rate** (subplot 1)
   - Analyzes impact of varying bit rates across different sampling frequencies
   - Shows measured vs theoretical BER curves

2. **BER vs Amplitude** (subplot 2)
   - Evaluates effect of signal amplitude on system performance
   - Dynamic SNR calculation based on amplitude variation

3. **BER vs Sampling Frequency** (subplot 3)
   - Studies effect of sampling frequency on demodulation performance
   - Tests across different signal amplitudes

4. **BER vs SNR** (subplot 4)
   - Primary performance analysis showing BER degradation with decreasing SNR
   - Validates theoretical predictions against simulation results

## Signal Visualization

Time-domain plots display:
- **Input Signal:** Generated random binary data
- **Modulated Signal:** FSK-modulated signal with carrier frequency (fc = 9 Hz)
- **Output Signal:** Demodulated and recovered data after noise corruption and demodulation

## Key Implementation Details

- **Demodulation:** Energy-based FSK demodulation with 0.5 threshold
- **RNG Seeds:** Multiple seeds (1, 6, 7) for reproducible yet varied test scenarios
- **Error Calculation:** `biterr()` function for precise BER computation
- **Noise Addition:** MATLAB `awgn()` function with 'measured' power specification

## Usage

Run the MATLAB script:
```matlab
WirelessProject.m
```

The simulation will:
1. Generate comprehensive console output with BER and PER metrics
2. Display performance comparison plots
3. Show time-domain signal representations
