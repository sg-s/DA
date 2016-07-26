
pHeader;



%% Reproducing the Naturalistic Stimulus
% In this document, I attempt to reproduce as closely as possible the naturalistic stimulus we use in fig 1. The problem with that stimulus is that it was generated using old MFCs with some arbitrarily chosen control signals, with bizarre valve arrangements, so we don't really understand what's going on, or how to modify it. Also, since it uses the Aalborg MFCs, it is less reproducible than what we can get using Alicats. 

%%
% This is what the stimulus that we want to reproduce looks like. In the following figure, we plot the stimulus (showing a small portion of it), together with regions where the stimulus in increasing. 

load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
S = mean(data(2).PID);
filt_window = 100;
dS = diff(filtfilt(ones(filt_window,1),filt_window,S));
dS(dS<0) = 0 ;
dS = filtfilt(ones(filt_window,1),filt_window,dS);
[ons,offs] = computeOnsOffs(dS>.5e-3);

rm_this = offs-ons<10;
ons(rm_this) = [];
offs(rm_this) = [];


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(S)
peak_S = NaN*ons;
for i = 1:length(ons)
	peak_S(i) = max(S(ons(i):offs(i)));
	plot(ons(i):offs(i),S(ons(i):offs(i)),'r')
end
xlabel('Time (10^4 samples /s)')
ylabel('PID (V)')
set(gca,'XLim',[2.3e5 3e5])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now that we have every putative valve on event, and the peak stimulus during that, we can export this data to an optimisation routine that attempts to match this. 




%% Version Info
%
pFooter;
