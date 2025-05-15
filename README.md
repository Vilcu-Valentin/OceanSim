# OceanSim
 
### Table of Content
- [Description](#description)
- [Documentation](#documentation)
- [Features](#features)
- [How to Run](#how-to-run)
- [Contributors](#contributors)

## Description

OceanSim is a real-time, physically based ocean surface for Unity.  
The water is generated from a Fast Fourier Transform (FFT) height field driven by a Phillips spectrum.
The simulation runs entirely on the CPU but is multi-threaded through Unity Jobs and Burst, and exposes utility functions for sampling height and normal so that boats or other rigid bodies can float convincingly on the dynamic surface.

## Documentation

## Features

- **Spectral Wave Synthesis** – Deep-water waves generated from a Phillips spectrum and evolved in the frequency domain  
- **2-D Inverse FFT Each Frame** – Real-time height, slope-x, and slope-z fields reconstructed on the CPU  
- **Unity Jobs + Burst** – Multi-core parallelism without requiring compute-shader support  
- **Tileable Ocean Plane** – Automatic instancing of surrounding tiles for an endless horizon  
- **Artist-Friendly Controls** – Wind speed, amplitude, chop, resolution, and colour gradients exposed in the inspector  

## How to Run

1. Clone the repository;

2. Open the project in Unity (recommended version: `2022.3.22f1`).

3. Open the `OutdoorsScene` in `Assets/Scenes`.

4. Press `Play` to start the simulation.

## Contributors

- [Valcu Valentin-Mihai](https://github.com/Vilcu-Valentin)
- [Iordache Alexandru-Stefan](https://github.com/IordacheAlexandruStefan)
- [Goidan Matei-Constantin](https://github.com/MateiGoidan)
- [Ilie Florin Alexandru](https://github.com/AlexFlorin21)
- [Ilie George Cristian](https://github.com/G3orge123)
