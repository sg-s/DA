
%% Pentyl Acetate Responses with ab2A
% In this document, we look at the responses of the ab2A neuron to pentyl acetate, which is pretty weird. 

pHeader;
dm = dataManager;

clearvars -except being_published dm
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('743b75637cedf714534e3e2e0814e776'),1);

cleanMSGdata

%% Data Overview
% First, I plot all the data, colour coded by the mean stimulus, to get a sense of what the data looks like:

c = parula(100);
mean_stim = mean(PID(a:z,:));
time = 1e-3*(1:length(PID));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:width(PID)
	ms = mean(PID(a:z,i));
	ci = 1+ceil(99*(ms - min(mean_stim))/max(mean_stim));
	plot(time,fA(:,i),'Color',c(ci,:))
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now, I plot only the responses from the highest mean stimuli, which show the strongest rebound after odor off:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:width(PID)
	if mean_stim(i) > .18
		ms = mean(PID(a:z,i));
		ci = 1+ceil(99*(ms - min(mean_stim))/max(mean_stim));
		plot(time,fA(:,i),'Color',c(ci,:))
	end
end
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% For these traces, we zoom in to see how the firing rate correlates with the stimulus. In the following plot, I have normalised the stimulus and the response so we can focus on the temporal structure of the two signals. The stimulus is shown in black and the firing rate in red.
do_these = (mean_stim > .18);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

		s = mean(PID(:,do_these),2); s = s - mean(s); s = s/std(s);
		r = mean(fA(:,do_these),2); r = r - mean(r); r = r/std(r);
		
		plot(time,s,'Color','k')
		plot(time,r,'Color','r')

xlabel('Time (s)')
set(gca,'XLim',[51 57])
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

% We clearly see that the firing rate massively increases when the odor turns off, and, more interestingly, the firing rate seems to be anti-correlated with the stimulus. 

%%
% When does this anti-correlation emerge? Obviously, at some point early in the game, the response is correlated with the stimulus. To look into this, I compute the cross-correlation of the stimulus and the firing rate in 2 second blocks and plot the extremal value of the cross-correlation as a function of time. In the following plot, positive extremal values indicate that increases in the stimulus correlate (with some lag) with increases in the response. Negative extremal values indicate that the increases in the stimulus correlate (with some lag) with DECREASES in the response. In other words, when the Y-values are greater than zero, the filter is positive, and when the Y-values are less than zero, the filter is negative. 

%%
% Also plotted are the location of the peak, which is the effective delay. 

% compute the cross-correlations
xc = NaN*fA;
d = NaN*fA;
bin_size = 2e3;
step_size = 50;


figure('outerposition',[0 0 1400 500],'PaperUnits','points','PaperSize',[1400 500]); hold on
subplot(1,3,1), hold on
clear l
for i = 1:width(PID)
	if do_these(i)
		for j = 2.5e3:50:57e3
			s = PID(j-bin_size:j+bin_size,i); s = s - mean(s); s = s/std(s);
			r = fA(j-bin_size:j+bin_size,i); r = r - mean(r); r = r/std(r);
			x = xcorr(r,s);
			x = x/(bin_size*2);

			xt = 1:length(x); xt = xt - mean(xt);

			if j == 6e3
				l(1) = plot(xt,x,'r');
			end
			if j == 50e3
				l(2) = plot(xt,x,'b');
			end


			if max(x) > -min(x)
				xc(j,i) = max(abs(x));
				[~,d(j,i)] = max(x); d(j,i) = d(j,i) - bin_size*2;
			else
				xc(j,i) = -max(abs(x));
				[~,d(j,i)] = min(x); d(j,i) = d(j,i) - bin_size*2;
			end
		end
	end
end

ylabel(['cross-correlation b/w ' char(10) 'stimulus and response'])
legend(l,{'t = 6s','t = 50s'})
title('Cross correlations')
xlabel('Lag (ms)')

subplot(1,3,2:3)
plot(time,xc,'k+')
xlabel('Time (s)')
ylabel(['Extremal value of cross-correlation' char(10) 'b/w stimulus and firing rate'])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end




%% Version Info
%

pFooter;

