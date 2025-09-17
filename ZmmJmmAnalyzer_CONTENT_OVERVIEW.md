# ZmmJmmAnalyzer - Content Overview

## Overview
The ZmmJmmAnalyzer is a comprehensive CMSSW analyzer package for studying Z → μμ and J/ψ → μμ decays. It contains tools for data analysis, submission scripts, and various plotting utilities.

## Directory Structure

### 1. `/miniAODmmmm/` - Main Analysis Package
This contains the core CMSSW analyzer for processing miniAOD data.

#### **Plugins** (`miniAODmmmm/plugins/`)
- **`miniAODmmmm.h`** - Header file for the main analyzer class
  - Defines the miniAODmmmm class inheriting from EDAnalyzer
  - Contains member variables for trigger handling, muon analysis, vertex fitting
  - Stores tree variables for Run, LumiBlock, Event information, muon kinematics

- **`miniAODmmmm.cc`** - Main analyzer implementation (473 lines)
  - Author: Nimmitha Karunarathna
  - Created: May 21, 2024
  - Processes PAT muons from miniAOD
  - Performs kinematic vertex fitting
  - Handles trigger information and Monte Carlo truth matching
  - Outputs ROOT trees with analysis variables

- **`GenStudyAnalyzer.h/.cc`** - Generator-level analysis
  - Similar structure for studying generator-level particles
  - Used for Monte Carlo studies and truth matching

- **`BuildFile.xml`** - CMSSW build configuration
  - Dependencies on Framework, PatCandidates, vertex fitting tools
  - ROOT and CLHEP integration

#### **Python Configuration** (`miniAODmmmm/python/`)
- **`miniAODmuonsRootupler_1.py`** - Main configuration script
  - Sets up GlobalTag for different data-taking periods (2022, 2023, 2024)
  - Currently configured for `141X_dataRun3_Prompt_v3`
  - Loads necessary CMSSW modules and services

- **`runGenStudyAnalyzer.py`** - Configuration for generator studies

#### **Tests** (`miniAODmmmm/test/`)
- **`test_catch2_miniAODmmmm.cc`** - Unit tests using Catch2 framework
- **`test_catch2_main.cc`** - Test runner
- **`BuildFile.xml`** - Test build configuration

### 2. `/dataSubmit/` - Data Processing and Submission Scripts

#### **Monte Carlo Studies**
- **`JPsiMuMu_MC/`** - J/ψ → μμ Monte Carlo analysis
  - Contains reconstructed events passing GEN conditions
  - Includes truth matching and muon soft IDs

#### **Data Production Versions**
- **`production_reduced_size_v1/` through `v6/`** - Evolution of data processing
- **`production_reduced_size_v6/`** - Latest version
  - Includes J/ψ rapidity×1000 for ATLAS-style analysis
  - Counts primary vertices
  - Flags true J/ψ vs combinatorial background

#### **Data Samples**
- **`ParkingDoubleMuonLowMass/`** - Low-mass dimuon parking data
- **`Muon/`** - Standard muon datasets
- **`ZeroBias/`** - Zero bias data for efficiency studies

#### **Timestamp Versions**
- **`with_timestamps_v1/` through `v3/`** - Time-stamped data processing

#### **Utility Scripts**
- **`checkFiles.sh`** - File integrity checking
- **`execute_crab_command.sh`** - CRAB job execution automation

### 3. `/scripts/` - Analysis and Plotting Tools

#### **Plotting Tools**
- **`plot_ratios/`** - Ratio plotting utilities
  - `plot_ratios.py` - Python script for creating ratio plots
  - `plot_ratios.ipynb` - Jupyter notebook interface
  - `plot_pt_eta.ipynb` - pT and η distribution plots
  - `utils/` - Supporting utility functions

- **`plot_eta/`** - Pseudorapidity distribution analysis
- **`Plot_ratio_fill/`** - Advanced ratio plotting with fills

#### **Analysis Scripts**
- **`cutBasedAnalyzer/`** - Cut-based analysis tools
- **`other/`** - Miscellaneous analysis scripts

#### **Luminosity Calibration**
- **`lumi_calibration_all_in_one/`** - Comprehensive luminosity tools
  - `scaling_to_brilcalc/` - Interface to BRIL calculator
  - Luminosity normalization and calibration utilities

## Key Features

### Analysis Capabilities
1. **Dimuon Reconstruction** - Z → μμ and J/ψ → μμ analysis
2. **Vertex Fitting** - Kinematic vertex reconstruction using KalmanVertexFitter
3. **Trigger Analysis** - HLT trigger path evaluation and prescale handling
4. **Monte Carlo Truth Matching** - Generator-level particle matching
5. **Luminosity Analysis** - Detailed luminosity measurement and calibration

### Data Processing
1. **Multi-year Support** - Handles 2022, 2023, 2024+ data
2. **Various Data Streams** - Parking, standard muon, zero bias datasets
3. **CRAB Integration** - Automated grid job submission
4. **Quality Assurance** - File checking and validation tools

### Physics Analysis
1. **Mass Spectrum Analysis** - J/ψ and Z boson mass reconstruction
2. **Kinematic Distributions** - pT, η, φ distributions
3. **Efficiency Studies** - Trigger and reconstruction efficiencies
4. **Background Studies** - Combinatorial background evaluation
5. **Cross-section Measurements** - Luminosity-normalized yields

## Usage
The analyzer is designed to be run within the CMSSW framework, processing miniAOD files to extract physics quantities of interest for muon pair studies. The output ROOT trees can be analyzed using the provided Python scripts and Jupyter notebooks.

## Recent Updates (as of 2024)
- Enhanced J/ψ analysis with rapidity cuts matching ATLAS analysis
- Improved background rejection with true J/ψ flagging
- Updated GlobalTag for 2024 data processing
- Comprehensive testing framework with Catch2 integration