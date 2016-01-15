# Fast Gain Control for Naturalistic Odor Detection in *Drosophila* 

This repository contains all code written for this project. Everything is written in MATLAB.

## Project Overview

### Code 

All the code here is written in MATLAB. 


| folder        | Description |
| -------       | ----------- |
| ./src 	    | contains core source files for some low-level calculations. Add to your path | 
| ./experiments | contains code for running the experiments. Needs kontroller (see below) |
| ./obsolete    | obsolete code, should be ignored unless you know what you're doing | 
| ./scripts 	| scripts to run various analyses. Run by calling them, or use makePDF (see below) |
| ./paper-figures | scripts that, when run, make figures for the paper. |


### Data 

You will need the data. The data is not part of this repository. You need to change paths everywhere in this repo. so that the data is referenced correctly. 

## Installation


### Download code

Install using my package manager:

```matlab
% copy and past this into your MATLAB shell
urlwrite('http://srinivas.gs/install.m','install.m')
install kontroller srinivas.gs_mtools
install spikesort da
install fitFilter2Data fitModel2Data
```

### install tag

Scripts assume you are running Mac OS X. Scripts also use `tag` to tag scripts when they successfully make a PDF. Install `tag` using 

```bash
brew install tag
```

This is totally not needed, but will throw errors if you don't install this. 

### Latex

You should have latex (Tex Live 2013+) installed. Check with:

```bash
latex --version
```

### makePDF

All scripts in this repository are written assuming that will be called using makePDF.m (part of srinivas.gs_mtools). For a tutorial on why I'm doing this, and why it matters, see [this repo](https://github.com/sg-s/awesome-matlab-notebook).

## Tests

To check that everything is working as it should be, you can run a test script that will attempt to run all scripts in `scripts` and in `paper-figures`. If this runs with no errors, making PDFs (see next section) should be OK. 

## Making PDFs

You can generate PDFs from any script (either in `scripts` or in `paper-figures`) using:

```matlab
makePDF(script_name.m)
```
Make sure you read this [tutorial](https://github.com/sg-s/awesome-matlab-notebook) to familiarize yourself with how makePDF works. 



