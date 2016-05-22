
%% Weber's Law generally observed
% In this document we look at a number of odorants and receptors and check if Weber's Law is true in the LFP and the firing rate. We do so by presenting Gaussian stimuli with increasing means. 


pHeader;
dm = dataManager;

%% 1-pentanol and ab3A
% In this section, we use 1-pentanol to stimulate the ab3A neuron. 
clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('29e79fc10c3dfd3c6cc9ad82a834f41b'),1);


cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- ab3A')

%% 1-pentanol and ab2A
% In this section, we use 1-pentanol to stimulate the ab2A neuron. 
clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('cd6753c0e4cf02895cd5e2c5cb58aa1a'),1);

cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- ab2A')

%% 2-butanone and ab2A
% In this section, we use 2-butanone to stimulate the ab2A neuron. 
clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('3ea08ccfa892c6545d74bbdaaa6cbee1'),1);

cleanMSGdata
makeMSGplots
suptitle('2-butanone -- ab2A')

%% butanal and ab3A
% In this section, we use butanal to stimulate the ab2A neuron. 
clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('5215c83836f0e4afa3060dab51081743'),1);

cleanMSGdata
makeMSGplots
suptitle('butanal -- ab3A')

%% ethanol and ab3A
% In this section, we use ethanol to stimulate the ab2A neuron. 
clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('2fd8239d788ef142be24f19e85e70427'),1);

cleanMSGdata
makeMSGplots
suptitle('ethanol -- ab3A')



%% Version Info
%
pFooter;
