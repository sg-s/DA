# Dynamical Adaptation

This document outlines all files in this project. Portions of this outline that are still to be completed are highlighted in *italics*.
## 1. Filter Analysis
`Question: How can we best back out a linear filter from a dataset of a stimulus and a response?`

	script to run: Filter_Analysis.m

This document shows demonstrates the validtiy of the algorithms I use to reconstruct the linear filter from input-ouput data. Specifically, it has the following sections: 

* Overview of different regularisation methods.
* Synthetic Data with No regularisation
* Synthetic Data 1: effect of regluarisation (here the stimulus is white noise)
* Synthetic Data 2: effect of regularisation (here the stimulus is filtered white noise)
* Real Data 1: effect of regularisation (1-octen-3-ol flickering stimulus data)
* Real Data 2: effect of regularisation (fast odor flickering stimulus data)
* Real Data 3: effect of regularisation ("Natural stimulus" data)


## 2. Effect of Correlation Time on Filters
`Question: Does (and how) the correlation time of the stimulus affect the shape of the recostructed filter?`

	script to run: EffectOfCorrrelationTime.m


## 3. Gain Analysis of 1-octen-3-ol flickering stimulus

`Question: Does ab3A rapidly modulate gain in response to 1-octen-3-ol?`

	script to run: GainAnalysisExample.m

In this document, we take one particular dataset (Carlotta's measuremnet of response of ab3 to 1-octen-3-ol) and check if a linear gain analysis can tell us if the neuron is adapting gain on a fast time scale. This document has the following sections: 

* Overview of the Data
* Statistics of the Stimulus and the Response
* Linear Fit to the Data
* Adding a output nonlinearity post-hoc
* Linear Model Gain Analysis: Fixed History Length
* Linear Model Gain Analysis: Varying History Length
* LN Model Gain Analysis: Fixed History Length
* LN Model Gain Analysis: Varying History Length




## 3. Analysis of "naturalistic stimulus"

## 4. Analysis of all of Carlotta's data

For most of Carlotta's data, the ORN is compeltely silent for most of the time, and the normal methods of gain analysis won't work so well. here we try to use other methods. 

## 5. Spike Frequency Adaptation Models