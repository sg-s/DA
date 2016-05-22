
%% Pentyl Acetate Responses with ab2A
% In this document, we look at the responses of the ab2A neuron to pentyl acetate, which is pretty weird. 

pHeader;


clearvars -except being_published
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/Users/sigbhu/ephys/sort_Me/pentyl-acetate/ab2',1);

cleanMSGdata

%% Data Overview
% The following plots show the PID stimulus and the LFP responses, normalised to show how the LFP  lag the stimulus. I've flipped the LFP upside down, so increases in LFP plotted here should correspond to increases in stimulus.

show_these_trials = [12 16 19 ];

a = 35e3;
z = 37e3;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
for i = 1:length(show_these_trials)
	subplot(1,length(show_these_trials),i); hold on
	s = PID(a:z,show_these_trials(i));
	s = s - mean(s);
	s = s/std(s);
	plot(s,'k')
	r = -LFP(a:z,show_these_trials(i));
	r = r - mean(r);
	r = r/std(r);
	plot(r,'r')
end




%%
% We note that there are two very weird aspects to this dataset. first, the firing lags are absolutely massive (~1s), or 10x as much as we would expect. Second, we see a very strange increase in the firing rate AFTER the odor is turned off. We have no explanation for this, and we confirm that this is odor-specific. The same sensillum did not show this behaviour with 1-pentanol

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(fA)
xlabel('Time (ms)')
ylabel('Firing Rate (Hz)')

%% Version Info
%

pFooter;

