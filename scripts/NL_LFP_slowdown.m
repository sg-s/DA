
pHeader;

%% Can NL models reproduce apparent LFP slowdown?
% In this document, I attempt to first fit the LFP data using a NL model, and then perform the slowdown analysis on the LFP, to see if the LFP does indeed slowdown. The following figure shows the LFP and the best-fit NL model. 

% get the data
load('/local-data/DA-paper/data-for-paper/fig7/nat-stim-ab3/2016_02_09_CS_F2_ab3_5.mat')
LFP = data(2).voltage(:,1:10:end)';
LFP = mean(LFP,2);
LFP = LFP - mean(LFP(1:5e3));
PID = data(2).PID(:,1:10:end)';
PID = mean(PID,2);
PID = PID - min(PID(1:5e3));

% fit the NL model to it
clear data p
data.response = LFP;
data.stimulus = [PID LFP];
p.k_D = 0.0835;
p.n = 1;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-3*(1:length(LFP));
plot(time,data.response,'k')
[R,K] = NLModel(data.stimulus,p);
plot(time,R,'r')
legend({'Data',['NL model, r^2 =' oval(rsquare(data.response,R))]},'Location','southeast')
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[0 70])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Now we perform the analysis of whether the LFP slows down or not, with both the real LFP and the NL model generated LFP. It looks like somehow, the NL model can reproduce what we see as a slowdown in the LFP. This means this is probably not a real effect. In addition, I plot the stimulus autocorrelation time (measured in 1 s windows), where the stimulus autocorrelation time is defined as the lag when the stimulus autocorrelation function drops below 1/e. 


min_acceptable_corr = .5;
min_acceptable_lag = 2;


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),-LFP(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','k','nbins',19);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),-R(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','r','nbins',19);

% find auto correlation time
act = NaN*PID;
mean_x = NaN*PID;
for i = 5e3:50:(length(PID)-5e3)
	x = PID(i-5e2:i+5e2);
	a = autocorr(x,length(x)-1);
	act(i) = find(a<1/exp(1),1,'first');
	mean_x(i) = mean(x);
end
rm_this = isnan(mean_x) | isnan(act);
mean_x(rm_this) = []; act(rm_this) = [];
plotPieceWiseLinear(mean_x,act,'Color','g','nbins',19);

xlabel('\mu_{Stimulus} in preceding 1s (V)')
ylabel('Lag (ms)')

legend({'Data','NL model','Stimulus auto corr.'},'Location','northwest')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Perhaps it's the stimulus correlations that cause this? What if we play back white noise into this model, and then repeat the analysis? I did this for three different NL models, that were identical except for their filters. I then plot the lag as a function of mean stimulus in the preceding 1s (exactly like before), and also overlay the location of the peak of the filter in the NL model. They line up nicely, as expected. 


clear p
p.k0 = 1;
p.tau = 100; % doesn't matter
p.B = 0; % no adaptation 

% parameters for filter 
p.tau1 = 10;
p.tau2 = 20;
p.n = 2;
p.A = 0;


% parameters for output scale
p.C = 1;


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

all_tau1 = [2 10 20];
for i = 1:length(all_tau1)
	% generate some new data
	S = randn(1e5,1);
	x = hill([1 .5 1],S);
	p.tau1 = all_tau1(i);
	K = filter_gamma2([],p);
	R = filter(K,1,x);

	[lag, mean_x] = findLagAndMeanInWindow(S(5e3:end-5e3),R(5e3:end-5e3),1e3,50);

	rm_this = isnan(lag) | isnan(mean_x);
	plotPieceWiseLinear(vectorise(mean_x(~rm_this)),vectorise(lag(~rm_this)),'Color','k','nbins',19);

	[~,loc] = max(K);
	plot([-5 5],[loc loc],'k--')
end


xlabel('Mean Stimulus (a.u.)')
ylabel('Lag (ms)')
title('NL model + white noise')
set(gca,'XLim',[-.15 .16])

set(gca,'YLim',[0 50])



prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% 
% What if I fit a NLN model to the firing rate data, and then repeat the analysis? In the following figure, I use a NLN model that I fit to the firing rate data, and treat that as data and repeat the analysis of the lag as a function of the mean stimulus. Both a NLN model and a NL model fit to the firing rate have a constant lag that does not depend on the mean stimulus, in contrast to the a NL model fit to the LFP data (see previous section). 

load('/local-data/DA-paper/data-for-paper/fig7/nat-stim-ab3/combined_data.ORNData','-mat')

clear data p
data.response = mean(od(3).firing_rate,2);
data.stimulus = [mean(od(3).stimulus,2) mean(od(3).firing_rate,2)];
data.stimulus(:,1) = data.stimulus(:,1) -min(min(od(3).stimulus(1:5e3,:)));
p.k_D = 0.0797;
p.n = 1.7561;
S = data.stimulus(:,1);
[R_NLN,K] = NLNmodel(data.stimulus,p);
filtertime = 1e-3*(1:length(K)) - .05;


% also generate NL model responses
x = hill([1 p.k_D p.n],S);
time = 1e-3*(1:length(S));
R_NL = convolve(time,x,K,filtertime);
R = data.response;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
[lag, mean_x, max_corr] = findLagAndMeanInWindow(S(5e3:end-5e3),R(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','k','nbins',19);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(S(5e3:end-5e3),R_NLN(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','r','nbins',19);
xlabel('\mu_{Stimulus} (V)')
ylabel('Lag (ms)')
legend({'Data','NLN model'})
set(gca,'YLim',[0 140])


subplot(1,2,2); hold on
[lag, mean_x, max_corr] = findLagAndMeanInWindow(S(5e3:end-5e3),R(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','k','nbins',19);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(S(5e3:end-5e3),R_NL(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','r','nbins',19);
xlabel('\mu_{Stimulus} (V)')
ylabel('Lag (ms)')
legend({'Data','NL model'})
set(gca,'YLim',[0 140])

prettyFig();
 
if being_published
	snapnow
	delete(gcf)
end


%% Showing slowdown in LFP and invariance of firing kinetics using Gaussian data
% What does this mean for the validity of this result? Cao et al. and Nagel et al. have showed that the LFP does slow down with adaptation. Martelli et al. have also shown that the firing rate does not. Can we still reconcile these two? We have another dataset where we have both LFP and firing rate measurement (the Gaussian stimulus dataset). In the following figure, I estimate lags using cross-correlations from this dataset, and compare lags of LFP and firing rate vs. the stimulus. 

%%
% In the following figure, I first estimate lags from the stimulus to the firing rate and the LFP using cross correlations, and plot that as a function of the mean stimulus (a). The error bars are standard error of the mean. Because we've been so screwed by cross correlations, I want to verify that I can actually see the LFP slowdown in the raw data, with no data analysis whatsoever. To check this, I first plot the LFP responses to odor offsets, (b), and then normalize these responses (c). In both cases, you can clearly see that the LFP slows down with increasing stimulus background (brighter colours). 

% get the gaussian data

% get the filter from the Gaussian stimuli 
clearvars -except being_published 
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);
v2struct(MSGdata)

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on

LFP_lags = NaN*paradigm;
LFP_max_corr = NaN*paradigm;
firing_lags = NaN*paradigm;
firing_max_corr = NaN*paradigm;
a = 35e3+1; z = 55e3;
chunk_size = 1e3;

% compute LFP and firing rate lags
for i = 1:length(paradigm)
	S = PID(a:z,i);
	X = -raw_LFP(a:z,i);
	R = fA(a:z,i);

	% reshape into chunks
	S = reshape(S,chunk_size,length(S)/chunk_size);
	X = reshape(X,chunk_size,length(X)/chunk_size);
	R = reshape(R,chunk_size,length(R)/chunk_size);

	S = bsxfun(@minus, S, mean(S));
	X = bsxfun(@minus, X, mean(X));
	R = bsxfun(@minus, R, mean(R));

	X_lag = NaN(chunk_size*2-1,size(S,2));
	R_lag = NaN(chunk_size*2-1,size(S,2));
	clear X_lag R_lag
	for j = 1:size(S,2)
		X_lag(:,j) = xcorr(X(:,j)/std(X(:,j)),S(:,j)/std(S(:,j)));
		R_lag(:,j) = xcorr(R(:,j)/std(R(:,j)),S(:,j)/std(S(:,j)));
	end
	X_lag = X_lag/chunk_size;
	X_lag = mean(X_lag,2);
	[LFP_max_corr(i),loc] = max(X_lag);
	LFP_lags(i) = loc - 1e3;

	R_lag = R_lag/chunk_size;
	R_lag = mean(R_lag,2);
	[firing_max_corr(i),loc] = max(R_lag);
	firing_lags(i) = loc - 1e3;

end

firing_lags(firing_lags<-100) = NaN; % obviously wrong

mean_stim = mean(PID(a:z,:));
for i = 1:max(paradigm)
	x = nonnans(firing_lags(paradigm==i));
	errorbar(mean(mean_stim(paradigm==i)),mean(x),sem(x),'k')
	x = nonnans(LFP_lags(paradigm==i));
	errorbar(mean(mean_stim(paradigm==i)),mean(x),sem(x),'r')
end
xlabel('\mu_{Stimulus} (V)')
ylabel('Lag (ms)')
legend({'Firing rate','LFP'},'Location','northwest')
set(gca,'XLim',[0 1.7],'YLim',[0 200])

subplot(1,3,2); hold on
xlabel('Time (s)')
title('LFP offset responses')
ylabel('\DeltaLFP (mV)')

subplot(1,3,3); hold on
title('LFP responses to odor offset')
xlabel('Time (s)')
ylabel('LFP change (norm)')

c = parula(11);
za = 54e3; zz = 57e3;
time = 1e-3*(1:length(PID));

for i = 1:max(paradigm)
	X = LFP(:,paradigm==i);
	for j = 1:width(X)
		X(:,j) = X(:,j) - mean(X(za:55e3,j));
	end
	X = nanmean(X,2);

	subplot(1,3,2); hold on
	plot(time(za:zz),X(za:zz,:),'Color',c(i,:))

	subplot(1,3,3); hold on
	X = X/max(X(za:zz));
	plot(time(za:zz),X(za:zz,:),'Color',c(i,:))

end

prettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;



