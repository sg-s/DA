# compiled feedback and criticism of carlotta's paper (we should avoid this)

1. Katherine Nagel thought that "slow" odors are only a consequence of paper dilution, and wonders if you still have slow odors if you do gas phase dilution.
2. Carlotta's paper contradicted Nagel & Wilson 
3. odor delivery system poorly described 
4. results shown only for one ORN and 2 odors
5. fast flickering stimulus pre-empts slow temporal patterning
6. claims of "similar" and "different" without statistical quantification 
7. how do Carlotta et al. see 6ms changes in filters if their window size of calculating firing rates is 100ms? 

# raw feedback that carlotta's paper got


## Neuron Reviewer 1

The authors' main claims are
i) different odors display different temporal dynamics 
ii) those dynamics are captured in the ORN responses and 
iii) that this could be a relevant feature in identifying different odors, or tracking plumes.  

They do an excellent job persuading me of the first two points, but the third one is more of an inference from their data, with little direct experimental evidence to support it.  Specifically, the authors do not show us whether flies can actually distinguish between 'fast' and 'slow' ORN activation, and whether this is related to the fly's ability to identify odors. (Aside from a role in specifying odor identity, the authors could have investigated whether flies' odor tracking behavior is dependent on whether they were following an odor with 'fast' vs 'slow' kinetics.  This is the most likely effect of different odor dynamics, and would be interesting (although perhaps not too surprising) to show.)  

Providing evidence for this third point - a role in odor identification - is the crucial aspect of this paper, since the first two points are interesting to those studying olfaction but do not reveal much new about the biology.  

The only evidence the authors present that odor dynamics contains information about odor identity is that different odors can be classified by the basic shape of their odor response (Fig 3c).  This shows that the dynamics are reliably different from one another, but does not establish whether flies use this information for odor identification purposes.  

If the authors want to claim that stimulus dynamics are important for odor identity, they need to address a few points.  The first is that different odors could likely exhibit exactly the same dynamics depending on  the airflow and concentration.  For example, a 'fast' odor in a slow total airflow might show kinetics similar to a 'slow' odor in fast airflow.  Secondly, the same odor in different airflows might show very different kinetics, but it seems likely that the odor is still identical to the flies.  Consequently, I find it difficult to believe that odor dynamics are important for the animal to identify monomolecular odors.  Isn't it more likely that the fly is sensitive to the dynamics in order to track the odor, but the dynamics don't influence the perception of odor identity?  In any case, some type of behavioral experiments would help us interpret these results as either relevant to tracking odor kinetics and/or odor identity. 

Another way of stating this point is that a corollary of their claim that stimulus dynamics impact odor identification is that identical odors delivered with distinct dynamics will be perceived as having a different identity, which seems implausible.  

Stimulus dynamics seem more likely to be important for identification of odor mixtures, as these authors previously reported (Su et al 2011). The current work is very well done, the data presented very clearly and this will be an important paper for the olfaction field, however I am not convinced that the dynamics of a monomolecular odor is part of what signals its identity, and don't see the results presented here to have more impact than the Su et al paper, unless they were to show more directly that flies use this information to identify or track an odor.  

Specific points:

page6 para2 Hierarchical clustering.  Needs to be described more thoroughly, I can't tell what was the input to the clustering algorithm, and particularly how the format of the input captured the dynamics of the responses.  

page6 para2 line6 'Clusters of ORN response can also be observed when performing PCA on the time-dependent activity (Fig 3a)'.  There is no plot showing any clusters, Fig3a is just the raw data of the odor dynamics.  Perhaps this was meant to be Supp Fig 3a?  Please clarify.

page7 para3 - Stimulus dynamics in the absence of airflow. The authors say they are mimicking the situation inside a flower, but surely the odor concentrations in the natural situation would be at steady state, rather than dynamic.  Here the dynamics arises because they have just introduced the odor.  

page12 para1.  Authors say they find small differences in dynamics of different ORNs to the same odor.  To me the differences in response to 1-octen-3-ol in ab3A and ab3B (Fig3) appear to be very significant: the ab3B response is square while the ab3A response is peaky and adapting.



## Neuron Reviewer #2:

Martelli et al. describe an analysis of Drosophila primary olfactory neurons responding to various odorants and found odors with different physical properties induce different response dynamics in ORNs. These response dynamics are odor specific and are maintained over different stimulus intensities by adaptation mechanisms. The authors suggest that odor discrimination may rely on these odorant dynamics as well as odorant molecular structure.

This work reveals an interesting aspect of Drosophila ORN biology-that odorant-induced spike rates normalized by peak response are similar for many odorants (over non-saturating concentration ranges).  This suggests this scaling transformation could be used to maintain an odor image over various stimulus concentrations.  This is a nice description of the phenomenon. Different odorants arrive with different kinetics, which is not generally appreciated, and since ORN responses correlate to the presence of odorant as measured by PID, this suggests odorant dynamics could be a novel discrimination dimension. However, I have a real problem believing this based on this work. There is no demonstration here that these mechanisms are actually used for this purpose by the fly.  There is no demonstration that manipulating response dynamics alone affects discrimination.  This is certainly approachable in this system.

Minor issues:
The gray bars in Figures 1 and 3 should be a little darker. Same for the error bars in 2F.
Page 2 third paragraph, the DasGupta reference. That flies with only one functional class of ORN can discriminate may not be true-this was done before the IR family of odorant receptors was appreciated.
All NR values in Figure 6 are less than 1. Does this mean they are all within the noise range and thus not significant?


## Neuron Reviewer #3: 

In the submitted manuscript, Carlotta et. al. presented a series of experiments directed towards understanding of the role of temporal features of odor stimuli for neuronal coding. In a well-controlled experimental setup, they measured the responses of olfactory receptor neurons (ORNs) in Drosophila to temporally variable odor stimuli. The crucial feature, which distinguishes this work from many other papers about ORNs responses to odors, is that the authors carefully measured the temporal profiles of the odor stimulus using a photo ionization detector (PID). That sets a high standard in the field. I found very interesting the authors' observation that the firing rate response of the receptor neurons scales up with concentration. In some range of concentrations, responses normalized by the peak amplitude look very similar. To my knowledge, this has not been reported before. Another interesting observation is that the linear filter for the LN model describing the transformation from odor temporal profile to firing rate receptor response is almost independent of odor identity (aside from a small temporal shift). Therefore, for any odor, the temporal dynamics of the stimulus, and not odor identity, defines the dynamics of the receptor response.  Both of these observations are important. 

However I felt that authors were trying to say more than that, and I could not follow their logic. The authors demonstrated that due to specific properties of individual odorants, their temporal dynamic is different and thus can be 'read' by receptor neuron. So, looking on the spiking response of receptor neurons one can discriminate which odorant has been delivered from a given odor delivery system, independently of odorant concentration. How does this observation relate to natural olfactory problems?  In all three cases of odor delivery used in this study, the temporal dynamics carry information about specific odor. Is this a general feature? Some odorants evaporate slowly, some faster. In the conditions described in the paper, stimulus dynamics suffice to discriminate between them. But it is not clear if it is a general property? Would a 'slow' odorant in a strong wind with higher degree of turbulence be different from 'fast' odorant in a weak wind?  If authors are trying to say that there is enough information in the temporal dynamics to identify an odor, they need to systematically change conditions of odor delivery and demonstrate that odor can be identified independently of these conditions. In all possible space of odor delivery conditions, the three conditions of odor delivery used in this study seem randomly chosen. Thus it is not clear how authors can generalize their result to the arbitrary situation, and make a case that the temporal dynamics of receptor activity is important for odor identification.  

To my opinion, the observations reported in this paper do not justify publication in Neuron, especially because the authors postulate an unconvincing relationship between these properties and neuronal coding.

Minor points:
1. In multiple cases authors mixed the notions of discrimination between two stimuli with identification of an odorant.  Page 2, par.3, reference to DasGupta and Waddell: "combinatorial map is not always required to encode odor identity". The fact that flies can discriminate two odors using one receptor type does not mean they can 'identify' odors. I paid attention to this statement, because the authors misuse this concept throughout the manuscript.

2. If I understand it correctly, the authors used 100 ms sliding window with 10 ms step. This is a relatively large time interval, and it was surprising to see an effect of 6 ms filter shift. I would appreciate seeing some temporal analysis of these results. Are filters stable for different widths of the sliding window?

3. Page. 4 par. 3, the description of individual neuron variability (Fig.1e) is not clear.

4. Page 7, par.3. The experiment with 'stimulus dynamics in the absence of air flow' actually does entail air flow. First, the PID itself,  creates a strong suction from the odor source to PID. The MiniPID (from Aurora Scientific) in its lowest flow range sucks ~0.7 l/min (or 12 cc/min), which definitely cannot be considered as 'the absence of air flow'. Second, the authors said, 'In this situation evaporation from the liquid source and diffusion likely dominate stimulus dynamics". Diffusion may play a role only at spatial scales much less than mm. On all larger scales odors are carried not by diffusion but always by convection. It is almost impossible to create motionless air in lab conditions, and it never happens in nature. 
Just because it difficult to measure the air movement, does not mean that it plays no role in stimulus dynamics.

## J. Neuro Reviewer 1


1. The authors results are at odds with another recent study (by Nagel and Wilson). The other study claims that differences
in response kinetics is largely due to receptor-ligand kinetics. The authors will have to explain their results in context of the other paper. Why are their results different? Is it differences in the kinetics of odor delivery? Is it the choice of odors? Have the authors tried the same odor-receptor pairs as the other paper? Without a clear explanation of the differences between the two studies, this becomes just another paper that muddles the field rather than providing clarity. My interpretation is that receptor-ligand kinetics also play a role in determining odor dynamics (see for example in Figure 3 1-octen-3-ol response kinetics is different in the two neurons).

2. It seems that the authors in this paper have achieved really fast odor kinetics. Or have they? Figure 1a - the kinetics of the raster and the PSTH don't match up at all. The spikes in the raster start with a 50 ms delay which is magically gone in the PSTH. Please explain this discrepancy. It is worth having a panel with the same figure expanded so that we can see from -0.5 s to 0.5 s clearly.

3. A better description of their odor delivery system and the performance of their odor delivery system is essential to understanding this study. A figure akin to Nagel and Wilson - supplementary figure 1 will be very helpful to understand what the authors are doing as well as to understand the differences between the two papers.

4. The authors have only demonstrated concentration independent dynamics in the "sensitivity range". But the title and abstract does not mention the "sensitivity range". The statements in the abstract and title are not supported by data. Also, based on figure 1b, the sensitivity range occurs only after about 70 Hz peak response. If one considers the amplification from ORN to PN, then most of the dynamic range of the PNs is devoted 
to ORN responses that are below the sensitivity range. I am convinced that the dynamics of different odors are very different but not that they are invariant. My view based on the author's results is that some odors under certain concentration have invariant dynamics. Also, the authors are overselling
the idea of odor dynamics as a signature of the odor. The spatial code changes so much with the concentration of odor anyway.
5. I think that the most important result from this paper is that odor dynamics are so dependent upon the identity of odors and that the differences are so large. But I don't understand why 1) dynamics are so dependent on the vapor pressure of the odor 2) why is tau_on and tau_off correlated and dependent on vapor pressure. The authors could perhaps explain this result to us.

Minor:
Pg4, 2nd paragraph - the authors mean Fig 5 instead of 4.
The authors may have forgotten to add the details of PCA and clustering performed in Figure 4.


## J. Neuro Reviewer 2

I think the idea of intrinsically "fast" and "slow" odorants is quite interesting and maps onto work by David Laing (not cited here) using identical terminology, to describe latency effects in olfactory psychophysics.  But for a variety of reasons I'm not sure this manuscript is ready for publication, including the fact that many of the ideas presented here were actually tested in prior publications (see comments below).  Most importantly, It seems to me that the conclusion, summarized in the abstract, "This suggests that a single response function can be associated with a single ORN and mediates the response to a large set of different odors", is both unjustified and wrong (or useless, if the "large set of different odors" isn't intended to mean nearly all odors).

Unjustified, because the small sample size makes it difficult to credit a broad conclusion.  Results from the model are based on the responses of only one neuron to only two odors that happen to elicit similar responses in that neuron, at least in the beginning of the response (Figure 3b - see methyl butyrate and octen-3-ol). Also, by turning the stimulus on and off at a very fast pace in Figure 6, the authors preempted the slow temporal patterning in the responses, and explicitly increased the dependence of the responses on stimulus dynamics, making it difficult to know how the responses would play out without the flickering presentation. 

Wrong, because the conclusions do not explain some published observations. Consider the responses of ORN1a in Figure1a of Raman et al, 2010, a paper cited by the authors, to HEX and GER. These responses start at nearly the same time, but for HEX, the ORN goes silent for more than 1 second, while for GER, it increases its firing rate for about the same period. Such large differences in the type of the response (excitatory versus inhibitory) cannot simply be explained by differences in the stimulus dynamics. Surely stimulus dynamics are important and affect the responses, but it seems there also are stimulus-timing-independent differences in temporal patterns of responses evoked by different odors within a neuron.

I would suggest to the authors to try to extend their prediction analysis to include multiple odors that induce a variety of responses (particularly in the beginning of the response) and that they vary the stimulus slow enough to allow slow temporal patterns to have an influence. Then, rather than claim that stimulus dynamics alone is needed to predict responses to different odors in a given neuron, they could try to study the relative properties of stimulus dynamics and odor-specific patterning. This approach was taken in Brown et al 2005, although for a somewhat different purpose. 

Minor concerns:

(1) pg 3: "We show that adaptation capabilities maintain ORN response dynamics independent of stimulus and background intensities...."  And pg 5:  "For stronger responses... the degree of adaptation was nearly constant."  The authors should cite Ito et al., 2009 - figure 5 shows how adaptation and saturation heavily compress the responses of moth ORNs.

(2) Figure 2 (described on pg 6): The significance of these results is difficult to interpret because they lack statistical correction for multiple comparison.

(3) pg 6: section titled "Stimulus dynamics depend on odor type..."  This point has already been demonstrated by Su et al (2011).  The present authors cite this paper but need to make clear what has already been demonstrated and what remains unknown. 

(4) pg 7: "Remarkably, ab3A responses to four odors ranging from "fast" to "slow" closely followed the concentration profiles measured by the PID (Fig 3a)."  I just don't see this resemblance.  The authors should find some way to quantify this important comparison.

(5) pg 7: "To see how easy it would be to discriminate between these odors on the basis of time-dependent activation..."  This exact idea is tested in Figure 1C of Raman et al., 2010), which should be cited here.

(6) pg 8 "In this situation evaporation from the liquid source and transport by convection likely dominate stimulus dynamics (there is also a small flow due to the PID suction."  The opposite sounds more likely to me, that stimulus dynamics would likely be dominated by suction from the PID.  Can the authors find a way to demonstrate this point rather than assuming it?

(7) pg 8: From their model  the authors claim that the response functions for two odors shown in Figures 6C (left) and 6D (left) are "different" but that the ones shown in Figures 6E (left) are "very similar". I don't see such a clear trend.

(8) pg 8: "This result suggests that physico-chemical properties of odor molecules can affect the correlation time of the odorant signals..."  What is meant by "correlation time"?

(9) pg 10: "These properties impose strong constraints on the adaptation mechanisms of ORNs..."   Not sure what this sentence means.

(10)  Figure 3: it appears overlapping panels sometimes cut off parts (x-axis) of other panels.



## Rachel Wilson:

1)	Several comments on an appropriate definition and qualification of the region where the response dynamics rescales: e.g. [This seems somewhat misleading to me. There are clearly intensity-dependent dynamics “within the response range” – see e.g. your Figure 3a, your Figure 5, and our Figure 6a-b. It’s only the bottom portion of the range, but it seems perverse to call this not “within the range”. I think it would be fairer to say “we find that normalized response dynamics are independent of stimulus intensity for a portion of the neuron’s dynamic range”, or even “for a large portion of the neuron’s dynamic range”.]


2)	However with few exceptions most previous studies assumed stimuli to be independent of odor type and considered only one or few concentrations of the stimulus. Rachel’s comment: “Not true – we took some pains in this area, and in Figure 8d of our paper we show explicitly a case where the measured dynamical properties of the stimulus were indistinguishable across odors, but there were nevertheless odor-dependent dynamics in the neural response.”


## Kathy Nagel


The first set of questions pertains to odor delivery.  You note that
you used two different delivery systems, one with a puff of air
generated by a valve upstream of the odor, and one with constant flow
through an odorized bottle and the valve located downstream of the
odor.  In the results you note that differences between odor stimuli
that were present with the puff configuration were sometimes absent in
the bottle configuration, and that ab3a showed silent periods after
odor stimulation during flicker (bottle?) but not puff.  Both of these
observations suggest that the bottle configuration may be faster than
the puff configuration.  Is this true in general?  Along related
lines, is there any way to speed up the delivery of a "slow" odor for
example by creating a saturated vapor and performing dilutions in air?
 This would seem to be a good test of whether the dynamics are truly
invariant once you control for stimulus dynamics.


We observed that a constant flow and downstream valve configuration,
with the valve located just before the main delivery tube, gave rise
to more reliable odor delivery dynamics.  Here are some PID traces
showing that our stimulus dynamics were independent of odor identity
for most odors we tested (including 2-butanone 1:1 0 ("2but1" and
isoamyl acetate 0.5 "isoAMp5").  The exception were cases where the
PID response was very small (e.g. fenchone 0.25 here).  Since the PID
response here was so small, and since the other odors had similar
dynamics, my conclusion was the that slower dynamics in this case
represented a limit on PID sensitivity rather than a true difference
in dynamics.  However, I could have been wrong about this.  (Sorry for
the V(mV) label this is something I grabbed from an old lab meeting)

We did observe that there was a slow decay in PID response amplitude
over trials that was more pronounced with some odors than with others.
 It seems likely that this could be due to differences in vapor
pressure as you describe, especially as it is more pronounced for
isoamyl acetate.


The second set of questions relate to what you call firing rate
adaptation— the fact that for many stimuli, the spike rate peaks and
then decays during the odor stimulus.  In our paper we tried to argue
that in general this "adaptation" was actually the spiking response to
two different phases of transduction— the onset (slope), and the
steady state.  Do you think this could account for why there is a
constant relationship between these two over a wide dynamic range?
(since the onset slope and the peak of the transduction current should
be related).


In some of your traces it appears as if the decay in the firing rate
becomes somewhat slower at higher concentrations than at lower
concentrations (eg fig 5a).  I observed this as well in some of my
data (unpublished) for three different concentrations of 2-butanone:

I never published this because I could never quite figure out what was
going on.  You can see that at the lowest concentration there is no
adaptation of the receptor potential (LFP) but that the spike rate
does show transient and sustained components corresponding to the
onset slow and steady-state portions.  The filter analysis suggests
that it may be the LFP-to-spike transition that changes with intensity
but there is also much more pronounced adaptation of the LFP at higher
concentrations.  I am curious how often you observe this type of
behavior and if you have any interpretation for it?

A third question is about the odor-dependent delays in the filters.
These seem puzzling to me—do they depend on the form or extent of the
regularization?

With regards to the manuscript there are I do not think our results
are really in strong contradiction as you also observe cases where an
identical odor stimulus can produce different dynamics in different
neurons (fig 5, 3rd row).  I think it is an excellent warning to the
field that it is important to measure your stimulus and an interesting
finding that often stimulus dynamics depend on vapor pressure.
However there are a couple of places where I feel like difference
between our papers are emphasized in ways that may be misleading.  For
example, the abstract does not acknowledge that the same odor stimulus
can produce different dynamics in different ORNs which is an important
limit to the conclusion that most ORN dynamics can be explained in
terms of the stimulus. 


Also in the discussion, where you note that we
observed a large range of dynamics from single ORNs in response to
different odorants I think it might be fair to note that differences
could arise from differences in odor delivery (placement of the valve,
constant flow, etc. for controlled odor delivery, proximity to
surfaces for plumes) rather than imputing that we did not measure our
stimulus (we did!)


There are several places in the manuscript where you note that the
same ORN filter can produce either "tonic" or "phasic" responses
depending on whether the delivered stimulus is fast or slow (p. 9 4th
full paragraph, p. 14, 2nd full paragraph). This seems analogous to
our argument that the same spike filter can produce either tonic or
phasic responses depending on the slope of the LFP (fig. 4 in our
paper).  The source of a different onset slope could be different
binding kinetics (as we suggest) or different stimulus dynamics (as
you show).  Perhaps it would be beneficial to the field to note the
similarities of these two arguments as well as their differences.

One more small point: the comparison of pb1 response dynamics between
your paper and our figure 6 seems misleading because you are looking
at spike rates and that figure shows LFPs.  As noted above (and as
argued in our paper) spike rates show dynamics that are distinct from
and have different concentration-dependence from LFPs.

Again thank you very much for sharing your manuscript with us and we
hope that our comments will be helpful.

