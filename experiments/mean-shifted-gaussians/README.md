# Mean Shifted Gaussians

Scripts and functions in this folder make stimuli that consist of mean shifted gaussians with similar variances. 

# How to make nice mean shifted Gaussians 

1. Make ControlParadigms using MakeMeanShiftedGaussians.m
2. configure Kontroller and the olfactometer as described in MakeMeanShiftedGaussians.m
3. run the experiments, gather stimuli using this ansatz control paradigm
4. load the data from (3), and run nonLinearWarpMSG on this. this script will save a new set of ControlParadigms with "warped" in the name that, when run, should give you nice distributions. 

# Tips and Tricks

* make sure you flush a lot of air through the odour vial immediately before the flicker. This eliminates headspace odourant, and gets you to steady state without going through transients, which may kill the neuron. 