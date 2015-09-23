# List of models being considered

## Linear-Nonlinear models

one filter block, passes through output nonlinearity. 

See:

* Chichilnisky, E. J. (2001). A simple white noise analysis of neuronal light responses. Network (Bristol, England), 12(2), 199–213.
* Dayan and Abbott (2001) Theoretical Neuroscience pg 48-51

## Nonlinear-Linear-Nonlinear models

See the paper on Hammerstein-Wiener models:

Hunter, I. W., & Korenberg, M. J. (n.d.). The identification of nonlinear biological systems: Wiener and Hammerstein cascade models. Biological Cybernetics, 55(2-3), 135–144. doi:10.1007/BF00341929

## Dynamical Adaptation Model

See:

* Clark, D. A., Benichou, R., Meister, M., & da Silveira, R. A. (2013). Dynamical Adaptation in Photoreceptors. PLoS Computational Biology, 9(11), e1003289. doi:10.1371/journal.pcbi.1003289

## Fast Gain Control Model 

My simplification of the DA model, where we directly fit n-th order corrections to a "gain-like" term that divisively contributes to the responses. 

## Spike-Frequency Adaptation Model 

See:

* Liu, Y.-H., & Wang, X.-J. (2001). Spike-frequency adaptation of a generalized leaky integrate-and-fire model neuron. Journal of Computational Neuroscience, 10(1), 25–45.

Here, they develop a point neuron model that has a buffering calcium variable. Show that this can account for slowing down of firing rate to steps of current. Also shows effects like forward masking, similar to predicton by DA model. 

## Amplification vs. Kinetics

See: 
* Arbib, M. A. (2004). The Handbook of Brain Theory and Neural Networks (pp. 94–97). chapter by Seung H. 

describes a one-ODE system that captures the relationship between timescale of response and the degree of amplification. This is of a very similar form to Liu and Wang. 



# Relevant Papers

## Vision

1. Baccus, S. A., & Meister, M. (2002). Fast and Slow Contrast Adaptation in Retinal Circuitry. Neuron. Identify two types of adaptation to contrasts in the retina in salamnders and rabbits. Fit LN models to their data. 
2. Rieke, F. (J Neuro). Temporal Contrast Adaptation in Salamander Bipolar Cells. Classic paper, studies the contrast adaptation of bipolar cells. Fit LN models. Looks at how responses changes when std. of signal changes, but mean remains same. 
3. Sakai, H. M., Ken-Ichi, N., & Korenberg, M. J. (1988). White-noise analysis in visual neuroscience. Visual Neuroscience, 1(03), 287. Shows how old these techniques are: LN models, guassian assumptions of data, and modelling photoreceptor, bipolar and retinal gainglion cells' responses. 
4. Kim, K. J., & Rieke, F. (n.d.). Temporal Contrast Adaptation in the Input and Output Signals of Salamander Retinal Ganglion Cells. Jneurosci.org. Another classic paper about contrast adaptation in retinal cells. LN model fits, mulitplie cascades, intermediate steps (as in Wilson) -- all here. 
5. Fairhall, A. L., Lewen, G. D., Bialek, W., & de Ruyter van Steveninck, R. R. (2001). Efficiency and ambiguity in an adaptive neural code. Nature. Nice paper where they consider how to design a neural code to a stimulus that is changing, and whose statistics are also evolving in time. Show that flies perform close to the theoretical maximum. 
6. Chander, D., & Chichilnisky, E. J. (n.d.). Adaptation to Temporal Contrast in Primate and Salamander Retina. Chichilnisky's paper on adaptation to temporal contrast in primate retinas. Shows that adaptation has different effects in different subtypes of retinal cells. 
7. DeWeese, M., & Zador, A. (1998). Asymmetric dynamics in optimal variance adaptation. Neural Computation, 10(5), 1179–1202. Interesting theoretical paper that looks at how to adapt to changing second-order statistics, e.g., change in input variance. make a prediction that optimal adaptation leads to asymmetric response dynamics. 

## Olfaction 

### Response of ORNs

1. Kaissling, K. E. (1998). Flux detectors versus concentration detectors: two types of chemoreceptors. Chemical Senses, 23(1), 99–111. Says there are two types of models for detectors (not really sure if they are) and has >20 parameters in his model that tries to predict ORN responses. No quantitative fits to data. 
2. Schuckel, J., Torkkeli, P. H., & French, A. S. (2009). Two Interacting Olfactory Transduction Mechanisms Have Linked Polarities and Dynamics in Drosophila melanogaster Antennal Basiconic Sensilla Neurons. Journal of Neurophysiology, 102(1), 214–223. doi:10.1152/jn.00162.2009. Record from ORNs, use a PID with a tracer gas, and look at the ORN's frequency response. 
3. Kim, A.J., Lazar, A.A. & Slutskiy, Y.B. System identification of Drosophila olfactory sensory neurons. J. Comput. Neurosci. published online doi:10.1007/s10827-010- 0265-0 (21 August 2010).
4. Geffen, M. N., Broome, B. M., Laurent, G., & Meister, M. (2009). Neural Encoding of Rapidly Fluctuating Odors - Geffen2009.pdf. Neuron, 61(4), 570–586. doi:10.1016/j.neuron.2009.01.021. Recrodings from PNs, but interesting nevertheless. LN models, consider 2D LN models to explain responses to rapidly fluctuating stimuli. 

### Adaptation and Dynamics of Response

1. Nagel, K. I., & Wilson, R. I. (2011). Biophysical mechanisms underlying olfactory receptor neuron dynamics. Nature Publishing Group, 14(2), 208–216. doi:10.1038/nn.2725. Their main result is that ORN responses can be decomposed into two steps: transduction and spike generation. They can model both, and do experiments where each step is modified individually. They study "adaptation" classically: using probe pulses on top of conditioning stimuli. 
2. Martelli, C., Carlson, J. R., & Emonet, T. (2013). Intensity invariant dynamics and odor-specific latencies in olfactory receptor neuron response. Journal of Neuroscience, 33(15), 6285–6297. doi:10.1523/JNEUROSCI.0426-12.2013

### Adaptation 
1. De Palo, G., Boccaccio, A., Miri, A., Menini, A., & Altafini, C. (2012). A Dynamical Feedback Model for Adaptation in the Olfactory Transduction Pathway. Biophysical Journal, 102(12), 2677–2686. doi:10.1016/j.bpj.2012.04.040. A model with dozens of parameters can explain  why the responses to a sequence of idenitcal pulses first decreases, then increases. Do experiments in newts, salamanders and mice. 
2. Reisert, J., & Matthews, H. R. (1999). Adaptation of the odour-induced response in frog olfactory receptor cells. The Journal of Physiology, 519(3), 801–813. doi:10.1111/j.1469-7793.1999.0801n. Show that frog ORNs respond to pulses (do a dose-response), then do that on top of a background. Show that there is some adaptation. 
3. Störtkuhl, K. F., Hovemann, B. T., & Carlson, J. R. (1999). Olfactory adaptation depends on the Trp Ca2+ channel in Drosophila. The Journal of Neuroscience. Show, using EAGs and beahvioural experiments, that Trp channels are necessary for adaptation to one odor. Time scale of adaptation here is >1min. 


### Responses to mixtures

1. Münch, D., Schmeichel, B., Silbering, A. F., & Galizia, C. G. (2013). Weaker Ligands Can Dominate an Odor Blend due to Syntopic Interactions. Chemical Senses. Show that there is signficant mixture processing at the periphery, and that ORN responses to a mixture is not a simple linear sum of the responses to the components. A "syntopic" model, proposed by Rospars, can explain their data. 
2. Rospars, J.-P. (2013). Interactions of Odorants with Olfactory Receptors and Other Preprocessing Mechanisms: How Complex and Difficult to Predict? Chemical Senses. This papers goes with the one above, where the "syntopic" model is explained. 
3. Rospars, J.-P., Lansky, P., Duchamp, A., & Duchamp-Viret, P. (2003). Relation between stimulus and response in frog olfactory receptor neurons in vivo. European Journal of Neuroscience. These guys fit a model with ~20 parameters to explain the responses of frog ORNs. The fits don't look very impressive. 