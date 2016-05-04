# Fast Gain Control for Naturalistic Odor Detection in *Drosophila* 

This repository contains all code written for this project. Everything is written in MATLAB.

## 1. Get the Code 

All the code is [here](https://github.com/sg-s/da) and [here](https://git.yale.edu/sg565/da) and is written in MATLAB. 


| folder        | Description |
| -------       | ----------- |
| ./src 	    | contains core source files for some low-level calculations. Add to your path | 
| ./experiments | contains code for running the experiments. Needs kontroller (see below) |
| ./obsolete    | obsolete code, should be ignored unless you know what you're doing | 
| ./scripts 	| scripts to run various analyses. Run by calling them, or use makePDF (see below) |
| ./paper-figures | scripts that, when run, make figures for the paper. |

Install using my package manager:

```matlab
% copy and past this into your MATLAB shell
urlwrite('http://srinivas.gs/install.m','install.m')
install kontroller srinivas.gs_mtools
install spikesort da
install fitFilter2Data fitModel2Data data-manager
```

## 2. Get the data 

You will need the data. The data is not part of this repository. You need to get the data and put it somewhere on your computer. This project uses [dataManager](https://github.com/sg-s/data-manager) to interface between code and the data. To use this, navigate to the folder where you have the data, and get dataManager to index the data:

```matlab
% you should be in the folder with the data
dm = dataManager;
dm.rehash;
```

## 3. Check your System

### install tag

Scripts assume you are running Mac OS X. Scripts also use `tag` to tag scripts when they successfully make a PDF. Install `tag` using 

```bash
brew install tag
```

This is totally not needed, but will throw errors if you don't install this. You'll have to remove all references to `tag` in the code.  

### Latex

You should have latex (Tex Live 2013+) installed. Check with:

```bash
latex --version
```

### makePDF

All scripts in this repository are written assuming that will be called using makePDF.m (part of srinivas.gs_mtools). For a tutorial on why I'm doing this, and why it matters, see [this repo](https://github.com/sg-s/awesome-matlab-notebook).


## 4. Making PDFs and figures

You can generate PDFs from any script (either in `scripts` or in `paper-figures`) using:

```matlab
makePDF(script_name.m)
```
Make sure you read this [tutorial](https://github.com/sg-s/awesome-matlab-notebook) to familiarize yourself with how makePDF works. 

You can generate individual PDF files for the figures using:

```matlab
makePDF(script_name.m,'--dirty')
```





