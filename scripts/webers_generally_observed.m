
%% Weber's Law generally observed
% In this document we look at a number of odorants and receptors and check if Weber's Law is true in the LFP and the firing rate. We do so by presenting Gaussian stimuli with increasing means. 


pHeader;


%% 1-pentanol and ab3A
% In this section, we stimulate the ab3A neuron with 1-pentanol.
clearvars -except being_published
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/1-pentanol/ab3',1);

cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- ab3A')

%% 1-pentanol and ab2A
% In this section, we stimulate the ab2A neuron with 1-pentanol.
clearvars -except being_published
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/1-pentanol/ab2',1);

cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- ab2A')


%% pentyl acetate and ab2A
% In this section, we stimulate the ab2A neuron with pentyl acetate
clearvars -except being_published
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/1-pentanol/ab3',1);

cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- ab3A')

%% 1-pentanol and ab2A
% In this section, we stimulate the ab2A neuron with 1-pentanol.
clearvars -except being_published
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/pentyl-acetate/ab2',1);

cleanMSGdata
makeMSGplots
suptitle('pentyl acetate -- ab2A')

%% EAG with 1-pentanol and 
% In this section, we measure the EAG by impaling the antenna with the electrode and stimulate with 1-pentanol.
clearvars -except being_published
[PID, LFP, ~, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/EAG/1-pentanol',1);

cleanMSGdata
makeMSGplots
suptitle('1-pentanol -- EAG')


%% EAG with ethyl acetate and 
% In this section, we measure the EAG by impaling the antenna with the electrode and stimulate with 1-pentanol.
clearvars -except being_published
[PID, LFP, ~, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/EAG/2ac',1);

cleanMSGdata
makeMSGplots
suptitle('ethyl acetate -- EAG')



%% Version Info
%
pFooter;
