swartz meeting day 1

# venky murthy 
figure ground segregation in the olfactory system

1. detection
2. discirimination
3. source localisation <--this is an interesting problem in high reynolds number envirnoments
4. learn association 

scene analysis and object recognition (in olfaction)

what is an olfactory object? some combination of elementary odourant cues?

specific exp question: how well can a mouse pick out an individual odourant in a mixture?

task: odour object recongition in variable backgrounds

need to recognise one target odourant, background odourant varies unpredictably 

--> the number of possible combinations of 1 target and 15 background odourants is very large (~16K)

the way the experiment is done:

when you train the mouse, you first present only the target, and then make the background more and more complex (where there are lots of compounds in the background)

humans are very poor in this task (<70% with 5 backgrounds)

result: performance degrades with increasing component # 
but why? 

idea: if we look at the glomerular activation pattern, maybe we can predict how and where failures occur. 

caveat: can't image all the glomeruli. can't image all the combinations. 

so they only look at single odourants. hmmm. they make a map of the glomeruli that activate for each odourant presented individually. 

assumption: mixture glomerular activation pattern is a linear sum of the individual glomerular representation pattern. 

they use this to construct a "masking index", and show that masking index correlates well with error rate in behaviour. 

## vivek jayaram

there are nice ring attractor models that can account for complex things like place and grin cells. but it's hard to test experimentally. 

insects are great navigators: bees, ants. 

see: neuser et al nature 2008

task: arena with one cool spot and hot spots everywhere. the cool spot is at a specific location marked with visual cues, and varies randomly. flies go to this spot, and learn that the visual markers sign where the cool spot is.

this is due to Central Complex (CX). ~3K neurons, sterotypic structure. a few sunapses from sensors and motors. 

they can pull out receptive fields of neurons in there using imaging and white noise stimulation. 

these cells are orientation selective. (only vertically oriented features). flies innately respond to vertical features. 

the region of activation in the ellipsoid body that activates in response to a stimulus in space varies from fly to fly -- but within a fly, there is a mappying from real space to position of activity in the ellipsoid body. 

this exactly matches where the fly is relative to the stripe

if there are two vertical stripes, the fly only keeps track of one stripe, and activity in the EB amtches it perfectly. but, if you look carefully, it looks like the bump of activity in the ellipsoid body, jumps from one stripe to the other. 

if you put it in a more complex visual scene, there is still only a single bump in activity. so it's not just a representation of the visual scene, but a representation fo what it cares about(?)

if there's no visual stimuli at all, then the activity in the Eb tracks the self-motion of the fly. 

this bump attractor works as working memory too -- standing flies remember where they were facing. 

# Yu Hu, Harvard University
Reduced descriptions of neural circuit dynamics 

motif analysis. motifs compared to random networks...but is this fair? can random networks be created by random processes?

questions:
- are small local motifs enough to predict network dynamics?
- which motifs are releveant to explain global dynamics? 

interest in that we can go from a filter to a network structure? 

## jonathan cannon

this is about homeotasis: syanptic or intrnisic 

how do you control neuronal activity by synaptic or intrinisc mechanisms?

always negative feedback. intrisic plasticity (which we call gain control) is a change in the tuning curve. he is talking about left-right shift in the tuning curve. 

why would you control both synaptic and intrinisc homeostasis (as opposed to just one)? 

we want to do this to optimise dynamic range of response (gain control)

so he's talking about multiplicative vs. additive gain control. uses mutual information to measure the quality of gain control. 

## julijana gjorjiena

crawling in drosophila larvae

homeostatic control of motor circuits (CPGs). activity dependent tuning indeveloping circuits. how do maps form in developing animals? how does topographic connectivity arise? you need neural activity to refine arboraisation. 

evidence: we see spontaneous activity in developing neural circuits. 

in larvae, there is a rhythmic contraction. but in embroys, there is no coordinated movement. there is a gradual development of muscle cooirdination as the embryo grows up. 

Q: is there a rhythmic motor pattern in adults, and is there a similar transition from uncoordinated movement to coordinated movement? 

## xiao-jing wang

flexible routing of information in the cerebral cortex

there is a hierarchy of timescales of response to stimulus -- V1 (and other early sensory areas) respond quickly, and as we go deeper into the brain, there is slower and slower responses. 

each cortical area gets inputs from thousands of other areas. you probably want to ignore most of them, and focus on what matters to the behaviour at hand

dendritic disinhibition as a mechanism for pathway-specific gating

## Bijan Pesaran 

flow of information in the posterior parietal cortex

what do correaltions between parts of the brain mean? strategy: manipualte bahviour, see how correlations between parts of the brain change. 

they can measure a single spike train and two LFPs in different parts of the brain in the monkey and they want to see how they're linked. 

can show that in a task where monkeys have to reach for and look at t target at the same trigger, the reach system slows down the saccade system. 

## Alex Vaughan

the representational stucture of the OFC

what is the natural way to encode the perceptual world? we go from stimulus space to some representiaton space. 

simniarly, what is the natural way to encode the congnitice world? how do we go from stimulus space to representations of decisions? 

you need to eoncode
p(choise was right)
E(reward size expectation)
E(value of choice)

ansatz: integrated value is the product of the confidence and the reward size 

see: spectral clustering to get local neighbourhoods. but they don't know how many clusters there are, and what the size of the nieghbhourhood is. to find this, they use a monte-carlo stability analysis. 

steve zucker says: that the eigenvalue of the laplacian should tell you the number of clusters. but this only works for fully isolated clusters. 

really nice analysis that is model free and found clusters in the data naturally, with no model prior. 

## Philip Sabes

State estimation for control

state estimation problem: how do we know where we are? so a combination of vision and proprioception. 

see: restricted boltzman machines



## Mate Lengyel







