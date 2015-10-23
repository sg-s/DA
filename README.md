# Dynamical Adaptation in ORNs

This repository contains all code written for this project. Everything is written in MATLAB.

## List of Figures

1. Gain control in ORN responses to naturalistic Inputs
2. White noise characterisation of gain 
3. Weber-like scaling 
4. Responses speed up with increasing mean stimulus
5. Fast gain control
6. Fast gain control is broadly observed
7. LFP also shows gain control 
8. Gain control with contrast changes
9. Gain change without odorant binding (ReaChR experiments)
10. DA model fits to data

### Supplementary Figures

1. No one nonlinearity can account for gain changes with mean stimulus
2. Fast gain control cannot be explained by a static output non-linearity 
3. 

## Installation

Install using my package manager:

```matlab
% copy and past this into your MATLAB shell
urlwrite('http://srinivas.gs/install.m','install.m')
install kontroller
install srinivas.gs_mtools
install spikesort
install da
install fitFilter2Data
install fitModel2Data
```

## tag

Scripts assume you are running Mac OS X. Scripts also use `tag` to tag scripts when they successfully make a PDF. Install `tag` using 

```bash
brew install tag
```

## Latex

You should have latex (Tex Live 2013+) installed. Check with:

```bash
latex --version
```

## Contents

* `experiments` contains scripts and functions to make experimental control paradigms. These control paradigms are config files for Kontroller
* `obp` contains some scripts used in analysis of OBP mutants, not related to the main project
* `paper-figures` contains scripts that will make all the figures in the paper. Each script operates on one dataset and may make more than one figure
* `janelia-poster` contains scripts to make some figures for a poster presented at the Swartz conference, together with the actual poster. 
* `scripts` contains all scripts that analyse various aspects of the problem. 
* `src` contains functions that are essential to this project, but not general-purpose enough to be moved to other projects. Make sure you add this to your path. 
* `pdfs` contains all PDFs made by any script in this repository

## Tests

To check that everything is working as it should be, you can run a test script that will attempt to run all scripts in `scripts` and in `paper-figures`. If this runs with no errors, making PDFs (see next section) should be OK. 

## Making PDFs

You can generate PDFs from any script (either in `scripts` or in `paper-figures`) using:

```matlab
makePDF(script_name.m)
```
