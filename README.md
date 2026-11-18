# FIR-IIR-Audio-Noise-Reduction
MATLAB-based design and comparative analysis of FIR and IIR digital filters for audio noise reduction, including Butterworth, notch, and window-based filtering techniques.
# FIR and IIR Audio Noise Reduction Using MATLAB

**Digital Signal Processing | MATLAB | FIR & IIR Filters | Audio Signal Processing**

## 📌 Project Overview

This project presents the design, implementation, and comparative analysis of **Finite Impulse Response (FIR)** and **Infinite Impulse Response (IIR)** digital filters for noise reduction in audio signals.

Developed as part of the **EE0857 Digital Signal Processing** course at **Istanbul Kültür University**, the project investigates how different digital filtering techniques can suppress unwanted noise while preserving the characteristics of the original music signal.

All filter designs and analyses are implemented using core MATLAB functions, without requiring the Signal Processing Toolbox.

## 🎯 Objectives

- Design FIR and IIR digital filters using mathematical principles.
- Remove 50 Hz power-line interference and high-frequency noise.
- Compare frequency responses, phase characteristics, and group delays.
- Evaluate filtering performance using Signal-to-Noise Ratio (SNR).
- Investigate computational complexity and filter efficiency.
- Implement fundamental DSP algorithms without specialized MATLAB toolboxes.

## ⚙️ Methodology

### 1. Signal Generation

A synthetic music signal is generated using eight musical notes from the C-major scale.

**Signal specifications:**

- Sampling frequency: 8,000 Hz
- Signal duration: 2.8 seconds
- Number of samples: 22,400
- Musical frequencies: 262–523 Hz, including harmonics

Three types of interference are introduced:

- 50 Hz power-line hum
- 3,000 Hz high-frequency sinusoidal interference
- Additive white Gaussian noise

### 2. IIR Filter Design

The IIR filtering system consists of:

**50 Hz Notch Filter**
- Second-order IIR filter
- Designed to suppress power-line interference

**Butterworth Low-Pass Filter**
- Sixth-order filter
- Cutoff frequency: 2,000 Hz
- Designed using analog prototype poles, frequency prewarping, and the bilinear transformation

Both causal and forward-backward filtering approaches are evaluated.

### 3. FIR Filter Design

The FIR filtering system uses the window method.

**High-Pass Filter**
- Blackman window
- 801 coefficients
- Designed to suppress low-frequency interference

**Low-Pass Filter**
- Hamming window
- 101 coefficients
- Designed to attenuate high-frequency interference

The cascaded FIR system contains 901 coefficients and exhibits linear-phase behavior with a constant group delay of 450 samples.

## 📊 Experimental Results

The following results were obtained from the synthetic audio experiment.

| Filtering Method | Output SNR (dB) | SNR Improvement (dB) |
|---|---:|---:|
| Noisy input | 2.27 | — |
| Causal IIR | 1.74 | -0.53 |
| Forward-backward IIR | 14.42 | +12.15 |
| Delay-compensated FIR | **14.92** | **+12.65** |

### Key Findings

- Both FIR and IIR designs successfully attenuated the targeted interference frequencies.
- FIR filtering achieved the highest measured SNR in the experiment.
- The FIR design preserved linear-phase characteristics at the expense of greater filter length and delay.
- IIR filtering required substantially fewer coefficients.
- Forward-backward IIR filtering improved the measured SNR but requires offline processing.
- The experiment demonstrates the trade-off between filtering performance, phase characteristics, and implementation complexity.

These findings are specific to the experimental signal and filtering conditions used in this project.

## 💻 MATLAB Implementation

The implementation includes:

- Discrete-time signal generation
- Convolution and LTI system analysis
- Discrete Fourier Transform (DFT)
- Sampling and aliasing demonstrations
- Z-transform and inverse Z-transform
- Pole-zero analysis
- Butterworth filter design
- Bilinear transformation
- FIR window-based filter design
- Frequency-response and group-delay analysis
- SNR evaluation
- Time-domain and frequency-domain visualization

### Requirements

- MATLAB R2016b or later
- No additional toolboxes required

### Running the Project

1. Download or clone this repository.
2. Open MATLAB.
3. Navigate to the project directory.
4. Open `dsp.m`.
5. Run the script.

The program generates numerical results, filter coefficients, and visualizations.

## 📚 Academic Background

The project is based on established digital signal processing principles, including:

1. J. G. Proakis and D. G. Manolakis, *Digital Signal Processing: Principles, Algorithms, and Applications*, 3rd ed., 1996.
2. A. V. Oppenheim, A. S. Willsky, and S. H. Nawab, *Signals and Systems*, 2nd ed., 1997.
3. E. Saatçı, *Digital Signal Processing, Lecture Notes 1–10*, Istanbul Kültür University, 2025.

## 🤖 Artificial Intelligence Disclosure

Claude (Anthropic) was used as an assistant during the development of this project.

Its contributions included:

- MATLAB code development and debugging
- Assistance with mathematical derivations
- Reimplementation of selected toolbox functionality
- Technical writing and English language editing
- Report organization and formatting

The project team reviewed the filter designs, mathematical derivations, and numerical results and verified the reported results by executing the MATLAB implementation.

## 👨‍💻 Authors

**Barış Erişen**  
**Ataberk Şensoy**

Department of Electrical and Electronics Engineering  
Faculty of Engineering  
Istanbul Kültür University  
Istanbul, Türkiye

**Academic Year:** 2026–2027

## 📄 Project Status

This repository contains an undergraduate Digital Signal Processing course project.

The work is intended for educational, technical demonstration, and further research purposes.

**This project has not been published as an IEEE conference or journal paper.**

## 🔬 Future Work

Potential extensions include:

- Testing with real-world audio recordings
- Evaluating different noise levels and signal types
- Comparing computational execution times
- Implementing filters on STM32 microcontrollers
- Investigating real-time embedded audio processing
- Expanding the experimental evaluation for potential academic publication
