# Kinematic Fit for Top Quark Mass Measurement

This project simulates the decay of a top quark ($t$) and utilizes a kinematic fit to determine its mass. The decay channel considered is $t \\rightarrow Wb \\rightarrow e\\nu b$. The simulation is performed for an hadronic collider environment, where only the transverse momenta of the charged particles (the b-quark and the electron) can be measured directly. The neutrino, which is an undetected particle, poses a challenge in reconstructing the event's kinematics.

## Project Overview

The core of this project is a C++ program that uses Monte Carlo methods and a kinematic fitting algorithm to reconstruct the top quark's mass. The primary goal is to demonstrate how to apply kinematic constraints to improve the precision of a mass measurement, even with missing information.

The program is structured in the following steps:

1.  **Event Generation**: Top quark decays are simulated using the `TGenPhaseSpace` class from the ROOT framework. The decay chain $t \\rightarrow Wb \\rightarrow e\\nu b$ is simulated, generating the four-momenta of the final state particles: a b-quark, an electron, and a neutrino.
2.  **Detector Simulation**: Realistic measurement uncertainties are applied to the generated particles' four-momenta. This step simulates the limited resolution of a particle detector. The uncertainties (resolutions) on momentum, azimuthal angle ($\\phi$), and polar angle ($\\theta$) are modeled using Gaussian distributions.
3.  Kinematic Fit: The central part of this project is the kinematic fit. This technique uses χ2 minimization to find the most probable four-momenta of the final state particles. It constrains the kinematic quantities of the neutrino, as it originates from the decay of a W boson. Since the W boson's mass is not represented by a simple Gaussian distribution, a Breit-Wigner term is added to the generalized chi-square function to perform the minimization correctly. The `TMinuit` class from the ROOT framework is used to perform the $\\chi^2$ minimization.
4.  **Data Analysis and Visualization**: Histograms are used to visualize the results, including the distribution of the reconstructed top quark mass before and after the kinematic fit.

## Code Structure

The program is written in C++ and uses the ROOT framework. Key components include:

  - `Particle` class: A custom class inheriting from `TLorentzVector` to store particle properties.
  - `fcn` function: The function to be minimized by `TMinuit`. It calculates the $\\chi^2$ value based on the measured and fitted parameters, including the kinematic constraints.
  - `ApplyResolution` function: Applies a Gaussian smearing to particle momenta and angles to simulate detector resolution.
  - `Event` function: Simulates the top quark decay and generates the true four-momenta of the final state particles.
  - `Kinematic_Fit` function: Sets up and runs the `TMinuit` minimization.
  - `main` function: The main loop that orchestrates the event generation, simulation, fitting, and data filling into histograms.

## Requirements

  - ROOT framework
  - A C++ compiler (e.g., g++)

## How to Compile and Run

A `makefile` is provided to simplify the compilation process.

1.  Make sure the ROOT environment is set up.
2.  Compile the code using the provided `makefile`:
    ```bash
    make
    ```
3.  Execute the program:
    ```bash
    ./main
    ```

The program will run the simulation and display several canvases with histograms of the results.
