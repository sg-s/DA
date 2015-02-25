# Dynamical Adaptation in ORNs

This repository contains all code written for this project. Everything is written in MATLAB. You might want to look at the following projects too, for additional code, etc:

* [Kontroller](https://github.com/sg-s/kontroller) I used this to acquire all data
* [spikesort](https://github.com/sg-s/spikesort) I used this to sort all spikes from the electro-physiological recordings
* [mtools](https://github.com/sg-s/srinivas.gs_mtools) Various useful functions in MATLAB, that all other repositories depend on.

## Installation

Install using my package manager:

```matlab
urlwrite('http://srinivas.gs/install.m','install.m')
install kontroller
install srinivas.gs_mtools
install spikesort
install da
```

## Make Paper figures

A script called `PaperFigures.m` makes all the figures for the paper. If you *don't* run Windows, and have [pdflatex](http://www.math.rug.nl/~trentelman/jacob/pdflatex/pdflatex.html) installed, you can make a PDF with all the figures using:

```
MakePDF('PaperFigures')
```

First runs of all scripts will be slow (~10 minutes) but subsequent runs will be much faster (~ seconds) thanks to [cache.m](https://github.com/sg-s/srinivas.gs_mtools/blob/master/docs/cache.md)

